use std::collections::HashSet;

use anyhow::Context;
use ndarray::{Array2, Axis, stack};
use ordered_float::OrderedFloat;
use rust_htslib::{self, bam::HeaderView};

pub mod plp_from_records;

pub static BASE2IDX: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'G' as usize] = 1;
    table[b'A' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b'C' as usize] = 4;

    table
};

pub struct PlpInfo {
    pub normed_count: Array2<f32>, // (4, step), 4: GATC
    pub major: Vec<usize>,
    pub minor: Vec<usize>,
}

impl PlpInfo {
    pub fn drop_low_ratio_ins_locus(self, ratio_thr: f32) -> Self {
        let low_ratio_locus = self
            .normed_count
            .axis_iter(Axis(1))
            .enumerate()
            .filter(|(_, ratios)| {
                ratios
                    .iter()
                    .max_by_key(|v| OrderedFloat(**v))
                    .copied()
                    .unwrap()
                    < ratio_thr
            })
            .map(|(idx, _)| idx)
            .collect::<HashSet<_>>();

        if low_ratio_locus.is_empty() {
            return self;
        }

        let normed_count = self
            .normed_count
            .axis_iter(Axis(1))
            .enumerate()
            .filter(|(idx, _)| !low_ratio_locus.contains(idx))
            .map(|(_, arr)| arr)
            .collect::<Vec<_>>();

        let normed_count = stack(Axis(1), &normed_count).unwrap();

        let major = self
            .major
            .into_iter()
            .enumerate()
            .filter(|(idx, _)| !low_ratio_locus.contains(idx))
            .map(|(_, v)| v)
            .collect::<Vec<_>>();
        let minor = self
            .minor
            .into_iter()
            .enumerate()
            .filter(|(idx, _)| !low_ratio_locus.contains(idx))
            .map(|(_, v)| v)
            .collect::<Vec<_>>();

        assert_eq!(major.len(), minor.len());
        assert_eq!(normed_count.shape()[1], major.len());
        Self {
            normed_count,
            major,
            minor,
        }
    }
}

pub struct BamHeaderSeqInfo {
    pub name: String,
    pub length: usize,
}

pub fn extract_seq_info_from_header(header_view: &HeaderView) -> anyhow::Result<BamHeaderSeqInfo> {
    let header = rust_htslib::bam::Header::from_template(header_view);
    let header_hashmap = header.to_hashmap();

    if !header_hashmap.contains_key("SQ") {
        anyhow::bail!("invalid bam header. SQ not found");
    }

    let target_seq_infos = header_hashmap.get("SQ").unwrap();
    if target_seq_infos.len() > 1 {
        anyhow::bail!("more than one target seqs. not supported now");
    }

    let first_target_seq = &target_seq_infos[0];

    // println!("first_target_seq: {:?}", first_target_seq);

    if !first_target_seq.contains_key("LN") {
        anyhow::bail!("invalid bam header. LN not found");
    }

    if !first_target_seq.contains_key("SN") {
        anyhow::bail!("invalid bam header. SN not found");
    }

    let ln_str = first_target_seq.get("LN").unwrap();

    let length = ln_str
        .parse::<usize>()
        .context(format!("parse {} to usize error", ln_str))?;

    let sn = first_target_seq.get("SN").unwrap();

    Ok(BamHeaderSeqInfo {
        name: sn.to_string(),
        length,
    })
}

pub fn pileup_counter(_bam_fname: &str) -> anyhow::Result<()> {
    // let mut index_reader = rust_htslib::bam::IndexedReader::from_path(bam_fname).unwrap();
    // index_reader.set_threads(num_cpus::get_physical()).unwrap();

    // let target_seq_info = extract_seq_info_from_header(index_reader.header())?;
    // let target_seq_length = target_seq_info.length;
    // // 6: A C G T GAP ins
    // let mut counter = Array2::<f32>::zeros((target_seq_length, 6));
    // let mut max_depth = 0;

    // for pileup in index_reader.pileup() {
    //     let pileup = pileup.unwrap();

    //     max_depth = max_depth.max(pileup.depth());
    //     let cur_position = pileup.pos() as usize;

    //     for alignment in pileup.alignments() {
    //         if alignment.is_refskip() {
    //             continue;
    //         }

    //         if alignment.is_del() {}
    //     }
    // }

    Ok(())
}

// #[cfg(test)]
// mod test {
//     use crate::pileup_counter::pileup_counter;

//     #[test]
//     fn test_plp_counter() {
//         let fpath = "/data-slow/qingke-deliver-data/20250818_240601Y0012_Run0003/Group_0/barcodes_reads_cons_gen_amplicon/Consensus/Bam/Group_0_Adaptor-barcode295-2.sort.bam";
//         pileup_counter(fpath);
//     }
// }
