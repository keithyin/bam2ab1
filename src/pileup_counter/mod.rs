use std::collections::HashSet;

use anyhow::Context;
use ndarray::{Array2, Axis, stack, s};
use ordered_float::OrderedFloat;
use rust_htslib::{self, bam::HeaderView};

use crate::utils::binary_search_lower_bound;

pub mod plp_from_records;

pub static BASE2IDX: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'G' as usize] = 1;
    table[b'A' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b'C' as usize] = 4;

    table
};

#[derive(Debug)]
pub struct PlpInfo {
    pub normed_count: Array2<f32>, // (4, step), 4: GATC . ratio
    pub major: Vec<usize>,
    pub minor: Vec<usize>,
}

impl PlpInfo {

    pub fn print_major(&self, col: usize) {
        let pos = binary_search_lower_bound(&self.major, &col);
        let start  = pos.saturating_sub(5);
        let end = (pos + 5).min(self.major.len());

        for cursor in start..end {
            println!("tt:{} -> {:?}", self.major[cursor], self.normed_count.slice(s![.., cursor]));
        }
        println!("-------------------------------------------------------------------------------")

    }

    // snp 5%, nh_indel: 10%, homo_indel: 45%
    pub fn modify_ratio(
        &mut self,
        seq: &[u8],
        snp_thr: f32,
        nh_indel_thr: f32,
        homo_indel_thr: f32,
    ) {
        self.normed_count
            .axis_iter_mut(Axis(1))
            .enumerate()
            .for_each(|(tt, mut gatc_ratio)| {
                let maj = self.major[tt];
                let mio = self.minor[tt];
                /* 因为是 gap left alignment。所以如果是比对的 ins region，major向后移动一下
                                     |
                                       |
                    consensus: A C G - T
                    smc1     : A C G G T
                */
                let maj = maj + mio.min(1);
                if maj >= seq.len() {
                    return;
                }

                let base = seq[maj];
                let base_idx = BASE2IDX[base as usize] as usize - 1;

                let mut is_homo = false;


                // if maj > 0 {
                //     is_homo |= seq[maj] == seq[maj - 1];
                // }
                if (maj + 1) < seq.len() {
                    let maj_base = seq[maj];
                    let cnt = ((maj+1)..(maj+4).min(seq.len())).into_iter().map(|pos| if seq[pos] == maj_base {1} else {0}).sum::<usize>();
                    is_homo = cnt == 3;
                }

                let gap_ratio = 1.0 - gatc_ratio.sum();

                let ins_region = self.minor[tt] > 0;

                if !is_homo {
                    // non homo

                    if !ins_region {
                        // non-homo non ins region

                        /*
                            non homo insertion
                                            |
                           consensus: A C G T
                           smc1     : A C G -
                        */

                        if gap_ratio < nh_indel_thr {
                            // gatc_ratio[base_idx] += gap_ratio;
                        }

                        for iter_idx in 0_usize..4 {
                            if iter_idx != base_idx {
                                if gatc_ratio[iter_idx] < snp_thr {
                                    gatc_ratio[base_idx] += gatc_ratio[iter_idx];
                                    gatc_ratio[iter_idx] = 0.0;
                                }
                            }
                        }
                    } else {
                        // non-homo ins region

                        /*
                            non homo insertion  ins region
                                            |
                           consensus: A C G - T
                           smc1     : A C G A T
                        */
                        if gap_ratio > (1.0 - nh_indel_thr) {
                            gatc_ratio.mapv_inplace(|_| 0.0);
                        }
                    }
                } else {
                    // homo

                    if !ins_region {
                        // homo non del region

                        /*
                            homo non insertion
                                          |
                           consensus: A C G G
                           smc1     : A C - G
                        */

                        if gap_ratio < homo_indel_thr {
                            // gatc_ratio[base_idx] += gap_ratio;
                        }

                        for iter_idx in 0_usize..4 {
                            if iter_idx != base_idx {
                                if gatc_ratio[iter_idx] < snp_thr {
                                    gatc_ratio[base_idx] += gatc_ratio[iter_idx];
                                    gatc_ratio[iter_idx] = 0.0;
                                }
                            }
                        }
                    } else {
                        // homo ins retion

                        /*
                            homo insertion  ins retion
                                            |
                                              |
                           consensus: A C G - T T
                           smc1     : A C G A T T
                        */
                        if gap_ratio > (1.0 - homo_indel_thr) {
                            gatc_ratio.mapv_inplace(|_| 0.0);
                        }
                    }
                }
            });
    }

    // 如果当前位点的 ACGT 的比例加起来<ratio_thr, 那么扔掉该位点
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

pub fn extract_seq_info_from_header(
    header_view: &HeaderView,
) -> anyhow::Result<Vec<BamHeaderSeqInfo>> {
    let header = rust_htslib::bam::Header::from_template(header_view);
    let header_hashmap = header.to_hashmap();

    if !header_hashmap.contains_key("SQ") {
        anyhow::bail!("invalid bam header. SQ not found");
    }

    let target_seq_infos = header_hashmap.get("SQ").unwrap();
    // if target_seq_infos.len() > 1 {
    //     anyhow::bail!("more than one target seqs. not supported now");
    // }

    let mut results = vec![];

    for seq_info in target_seq_infos {
        if !seq_info.contains_key("LN") {
            anyhow::bail!("invalid bam header. LN not found");
        }

        if !seq_info.contains_key("SN") {
            anyhow::bail!("invalid bam header. SN not found");
        }

        let ln_str = seq_info.get("LN").unwrap();

        let length = ln_str
            .parse::<usize>()
            .context(format!("parse {} to usize error", ln_str))?;

        let sn = seq_info.get("SN").unwrap();

        results.push(BamHeaderSeqInfo {
            name: sn.to_string(),
            length,
        });
    }

    Ok(results)
}

#[cfg(test)]
mod test {
    use crate::pileup_counter::PlpInfo;

    #[test]
    fn test_plp_info() {

        // let mut plp_info = PlpInfo{normed_count: }
    }
}
