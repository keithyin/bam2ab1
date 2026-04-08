use std::{collections::HashMap, io::Write, path};

use bam2ab1::{
    ab1::{
        data_process::{Plp2Ab1WithDeletionShrink, Plp2Ab1WithInsIdentifier},
        plp2ab1::TPlp2Ab1,
    },
    pileup_counter::{
        BamHeaderSeqInfo, extract_seq_info_from_header, plp_from_records::plp_from_records,
    },
};
use clap::Parser;
use gskits::{
    ds::ReadInfo,
    fastx_reader::{fasta_reader::FastaFileReader, read_fastx},
};
use rust_htslib::bam::{Read, Record};
use time;
use tracing;
use tracing_subscriber;

#[derive(Debug, Parser, Clone)]
#[command(version, about, long_about = None)]
struct Cli {
    #[arg(long = "bam", help = "input bam file (aligned)")]
    pub bam: String,

    #[arg(long = "ref", help = "reference fasta file")]
    pub reference: String,

    #[arg(short = 'o', help = "output ab1 file path")]
    pub output_fpath: Option<String>,

    #[arg(long = "width")]
    pub base_width: Option<usize>,

    #[arg(long = "chunkSize")]
    pub chunk_size: Option<usize>,

    #[arg(long = "ovlpSize")]
    pub ovlp_size: Option<usize>,

    #[arg(
        long = "insIdent",
        help = "insertion region identifier. only accept N or acgt"
    )]
    pub insertion_region_identifier: Option<String>,

    #[arg(long = "bamReadThreads", default_value_t = 1)]
    pub bam_reader_threads: usize,
}

impl Cli {
    pub fn get_output_filepath(&self) -> String {
        if let Some(ref fpath) = self.output_fpath {
            return fpath.clone();
        }

        let bam_path = std::path::Path::new(&self.bam);
        let root = bam_path.parent().unwrap_or(std::path::Path::new("."));
        let bam_stem = bam_path.file_stem().unwrap().to_str().unwrap();

        root.join(format!("{}.ab1", bam_stem))
            .to_str()
            .unwrap()
            .to_string()
    }
}

fn build_plp2ab1_transformer(cli: &Cli) -> Box<dyn TPlp2Ab1> {
    match cli.insertion_region_identifier.as_ref().map(|v| v.as_str()) {
        None => Box::new(Plp2Ab1WithDeletionShrink::new()),
        Some("N") => Box::new(Plp2Ab1WithInsIdentifier::new(
            bam2ab1::ab1::data_process::InsertionRegionIdentifier::N,
        )),

        Some("acgt") => Box::new(Plp2Ab1WithInsIdentifier::new(
            bam2ab1::ab1::data_process::InsertionRegionIdentifier::Acgt,
        )),
        Some(a) => {
            panic!("invalid param for insIdent '{}',  try N or acgt", a);
        }
    }
}

fn main() {
    let time_fmt = time::format_description::parse(
        "[year]-[month padding:zero]-[day padding:zero] [hour]:[minute]:[second]",
    )
    .unwrap();
    let time_offset =
        time::UtcOffset::current_local_offset().unwrap_or_else(|_| time::UtcOffset::UTC);
    let timer = tracing_subscriber::fmt::time::OffsetTime::new(time_offset, time_fmt);

    tracing_subscriber::fmt::fmt().with_timer(timer).init();

    let cli = Cli::parse();

    let output_fpath = cli.get_output_filepath();

    let output_path = path::Path::new(&output_fpath);
    let output_root = output_path.parent().unwrap();
    let output_stem = output_path
        .file_stem()
        .unwrap()
        .to_str()
        .unwrap()
        .to_string();

    let fasta_reader = FastaFileReader::new(cli.reference.clone());

    let fasta_records = read_fastx(fasta_reader);

    let fasta_records = fasta_records
        .into_iter()
        .map(|rec| (rec.name.clone(), rec))
        .collect::<HashMap<String, ReadInfo>>();

    if fasta_records.is_empty() {
        tracing::error!("empty fasta file. {}", cli.reference);
        return;
    }

    let mut bam_reader = rust_htslib::bam::Reader::from_path(&cli.bam).unwrap();
    bam_reader.set_threads(cli.bam_reader_threads).unwrap();

    let bam_header_sq_records = extract_seq_info_from_header(bam_reader.header())
        .unwrap()
        .into_iter()
        .map(|seq_info| (seq_info.name.clone(), seq_info))
        .collect::<HashMap<String, BamHeaderSeqInfo>>();

    let header = rust_htslib::bam::Header::from_template(bam_reader.header());
    let bam_header = rust_htslib::bam::HeaderView::from_header(&header);

    if fasta_records.len() != bam_header_sq_records.len() {
        tracing::error!(
            "fasta bam not match. The amount of query data in the FASTA file does not match the number of reference sequences in the BAM file."
        );
        panic!(
            "fasta bam not match. The amount of query data in the FASTA file does not match the number of reference sequences in the BAM file."
        );
    }

    fasta_records.iter().for_each(|(name, info)| {
        if let Some(bam_header_info) = bam_header_sq_records.get(name) {
            if info.seq.len() != bam_header_info.length {
                tracing::error!("fasta bam not match. seq_len not match");
                panic!("fasta bam not match");
            }
        } else {
            tracing::error!("fasta bam not match. query name not found in bam");
            panic!("fasta bam not match. query name not found in bam");
        }
    });

    let mut target_records: HashMap<String, Vec<Record>> = HashMap::new();

    let mut record = Record::new();
    loop {
        if let Some(Ok(_)) = bam_reader.read(&mut record) {
            let t_name =
                String::from_utf8(bam_header.tid2name(record.tid() as u32).to_vec()).unwrap();
            if target_records.contains_key(&t_name) {
                target_records.get_mut(&t_name).unwrap().push(record);
            } else {
                target_records.insert(t_name, vec![record]);
            }
            record = Record::new();
        } else {
            break;
        }
    }

    let plp2ab1_transformer = build_plp2ab1_transformer(&cli);

    target_records.into_iter().for_each(|(name, mut records)| {
        records = records
            .into_iter()
            .filter(|record| {
                !record.is_unmapped() && !record.is_secondary() && !record.is_supplementary()
            })
            .collect();

        if records.is_empty() {
            tracing::warn!("empty records",);
            return;
        }

        let bam_seq_info = bam_header_sq_records.get(&name);
        if bam_seq_info.is_none() {
            tracing::warn!("empty records",);
            return;
        }

        let seq_info = bam_seq_info.unwrap();
        let fasta_seq_info = fasta_records.get(&name);
        if fasta_seq_info.is_none() {
            tracing::warn!("empty records",);
            return;
        }
        let fasta_seq_info = fasta_seq_info.unwrap();

        let reference_sequence = &fasta_seq_info.seq;

        let chunk_size = cli.chunk_size.unwrap_or(reference_sequence.len());
        let ovlp_size = cli.ovlp_size.unwrap_or(0);

        let mut window_start = 0;
        while window_start < reference_sequence.len() {
            let window_end = window_start + chunk_size;
            let window_end = if (window_end + ovlp_size) > reference_sequence.len() {
                reference_sequence.len()
            } else {
                window_end
            };

            let mut plp_info = plp_from_records(&records, window_start, window_end);

            // println!("{:?}", plp_info.normed_count);

            // plp_info.print_major(3);

            plp_info.modify_ratio(
                &reference_sequence.as_bytes()[window_start..window_end],
                0.05,
                0.1,
                0.45,
            );

            // println!("{:?}", plp_info.normed_count);

            // panic!("");

            // plp_info.print_major(3);

            let plp_info = plp_info.drop_low_ratio_ins_locus(0.02);

            // println!("{:?}", plp_info.normed_count);
            // panic!();

            // plp_info.print_major(3);
            let peak_width = if cli.base_width.is_none() {
                Some(
                    ((u16::MAX) as usize / (window_end - window_start) - 1)
                        .min(20)
                        .max(3),
                )
            } else {
                cli.base_width
            };
            tracing::info!("peak_width: {peak_width:?}");

            let ab1_file = plp2ab1_transformer.transform(
                &plp_info,
                &reference_sequence[window_start..window_end],
                peak_width,
                Some(seq_info.name.clone()),
            );

            let idx = if name.contains("|") {
                name.split_once("|").unwrap().0.to_string()
            } else {
                name.clone()
            };

            let o_path = output_root.join(format!(
                "{}.{}.{}-{}.ab1",
                output_stem, idx, window_start, window_end
            ));
            match std::fs::File::create(&o_path) {
                Ok(mut file) => {
                    if let Err(err) = file.write_all(&ab1_file.to_bytes()) {
                        tracing::error!("write to file error. {}", err);
                    } else {
                        tracing::info!("Success, write to : {:?}", o_path);
                    }
                }
                Err(err) => {
                    tracing::error!("ERROR: create file error: {}", err);
                }
            }

            if window_end >= reference_sequence.len() {
                break;
            }
            window_start = window_end - ovlp_size;
        }
    });
}
