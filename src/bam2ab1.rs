use std::{collections::HashMap, io::Write, path};

use bam2ab1::{
    ab1::data_process::transform_plp_info_2_ab1_data_with_deletion_shrink,
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

#[derive(Debug, Parser, Clone)]
#[command(version, about, long_about = None)]
struct Cli {
    #[arg(long = "bam", help = "input bam file (aligned)")]
    pub bam: String,

    #[arg(long = "ref", help = "reference fasta file")]
    pub reference: String,

    #[arg(short = 'o', help = "output ab1 file path")]
    pub output_fpath: Option<String>,
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

fn main() {
    let cli = Cli::parse();

    let output_fpath = cli.get_output_filepath();

    let output_path = path::Path::new(&output_fpath);
    let output_root = output_path.parent().unwrap();
    let output_stem = output_path.file_stem().unwrap().to_str().unwrap().to_string();

    let fasta_reader = FastaFileReader::new(cli.reference.clone());

    let fasta_records = read_fastx(fasta_reader);
    let name2idx = fasta_records
        .iter()
        .enumerate()
        .map(|(idx, record)| (record.name.clone(), idx))
        .collect::<HashMap<_, _>>();

    let fasta_records = fasta_records
        .into_iter()
        .map(|rec| (rec.name.clone(), rec))
        .collect::<HashMap<String, ReadInfo>>();

    if fasta_records.is_empty() {
        eprintln!("ERROR: empty fasta file. {}", cli.reference);
        return;
    }

    let mut bam_reader = rust_htslib::bam::Reader::from_path(&cli.bam).unwrap();
    bam_reader.set_threads(num_cpus::get_physical()).unwrap();

    let bam_header_sq_records = extract_seq_info_from_header(bam_reader.header())
        .unwrap()
        .into_iter()
        .map(|seq_info| (seq_info.name.clone(), seq_info))
        .collect::<HashMap<String, BamHeaderSeqInfo>>();

    let header = rust_htslib::bam::Header::from_template(bam_reader.header());
    let bam_header = rust_htslib::bam::HeaderView::from_header(&header);

    assert_eq!(
        fasta_records.len(),
        bam_header_sq_records.len(),
        "fasta bam not match"
    );

    fasta_records.iter().for_each(|(name, info)| {
        if let Some(bam_header_info) = bam_header_sq_records.get(name) {
            assert_eq!(
                info.seq.len(),
                bam_header_info.length,
                "fasta bam not match"
            );
        } else {
            panic!("fasta bam not match");
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

    target_records.into_iter().for_each(|(name, mut records)| {
        records = records
            .into_iter()
            .filter(|record| {
                !record.is_unmapped() && !record.is_secondary() && !record.is_supplementary()
            })
            .collect();

        if records.is_empty() {
            eprintln!("WARN: empty records",);
            return;
        }

        let bam_seq_info = bam_header_sq_records.get(&name);
        if bam_seq_info.is_none() {
            eprintln!("WARN: empty records",);
            return;
        }

        let seq_info = bam_seq_info.unwrap();
        let fasta_seq_info = fasta_records.get(&name);
        if fasta_seq_info.is_none() {
            eprintln!("WARN: empty records",);
            return;
        }
        let fasta_seq_info = fasta_seq_info.unwrap();

        let reference_sequence = &fasta_seq_info.seq;

        let mut plp_info = plp_from_records(&records, seq_info.length);
        plp_info.modify_ratio(reference_sequence.as_bytes(), 0.05, 0.1, 0.45);
        let plp_info = plp_info.drop_low_ratio_ins_locus(0.02);

        let peak_width = Some(20);
        let ab1_file = transform_plp_info_2_ab1_data_with_deletion_shrink(
            &plp_info,
            reference_sequence,
            peak_width,
            Some(seq_info.name.clone()),
        );

        let idx = name2idx.get(&name);
        if idx.is_none() {
            eprintln!("WARN: empty records",);
            return;
        }

        let idx = idx.unwrap();


        let o_path = output_root.join(format!("{}.{}.ab1", output_stem, idx));
        match std::fs::File::create(&o_path) {
            Ok(mut file) => {
                if let Err(err) = file.write_all(&ab1_file.to_bytes()) {
                    eprintln!("ERROR: write to file error. {}", err);
                } else {
                    println!("Success, write to : {:?}", o_path);
                }
            }
            Err(err) => {
                eprintln!("ERROR: create file error: {}", err);
            }
        }
    });
}
