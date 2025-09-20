use std::io::Write;

use bam2ab1::{
    ab1::data_process::transform_plp_info_2_ab1_data,
    pileup_counter::{extract_seq_info_from_header, plp_from_records::plp_from_records},
};
use clap::Parser;
use gskits::fastx_reader::{fasta_reader::FastaFileReader, read_fastx};
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

    let fasta_reader = FastaFileReader::new(cli.reference.clone());
    let fasta_records = read_fastx(fasta_reader);
    assert!(
        fasta_records.len() == 1,
        "Only support one reference sequence"
    );

    let reference_sequence = &fasta_records[0].seq;

    let mut bam_reader = rust_htslib::bam::Reader::from_path(&cli.bam).unwrap();
    bam_reader.set_threads(num_cpus::get_physical()).unwrap();
    let seq_info = extract_seq_info_from_header(bam_reader.header()).unwrap();
    assert_eq!(
        reference_sequence.len(),
        seq_info.length,
        "Reference length mismatch"
    );

    assert_eq!(
        fasta_records[0].name, seq_info.name,
        "Reference name mismatch"
    );

    let mut record = Record::new();
    let mut records = vec![];
    loop {
        if let Some(Ok(_)) = bam_reader.read(&mut record) {
            records.push(record);
            record = Record::new();
        } else {
            break;
        }
    }

    records = records
        .into_iter()
        .filter(|record| {
            !record.is_unmapped() && !record.is_secondary() && !record.is_supplementary()
        })
        .collect();
    let plp_info = plp_from_records(&records, seq_info.length);
    let ab1_file =
        transform_plp_info_2_ab1_data(&plp_info, reference_sequence, Some(seq_info.name.clone()));

    let mut file = std::fs::File::create(&output_fpath).unwrap();
    file.write_all(&ab1_file.to_bytes()).unwrap();
    println!("write done: {}", output_fpath);
}
