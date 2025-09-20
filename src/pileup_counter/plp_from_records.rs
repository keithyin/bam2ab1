use std::collections::HashMap;

use gskits::gsbam::{
    bam_record_ext::BamRecordExt, plp_counts_from_records::compute_max_ins_of_each_ref_position,
};
use ndarray::{Array2, ArrayViewMut1, Axis, concatenate, s};
use rust_htslib::bam::{Record, ext::BamRecordExtensions};

use crate::pileup_counter::BASE2IDX;

pub fn plp_from_records(records: &Vec<Record>, target_len: usize) -> super::PlpInfo {
    let major_pos_ins =
        compute_max_ins_of_each_ref_position(&records, Some(0), Some(target_len), None);

    let mut major_pos_ins_vec = major_pos_ins
        .iter()
        .map(|(&k, &v)| (k as usize, v as usize))
        .collect::<Vec<_>>();
    major_pos_ins_vec.sort_by_key(|v| v.0);
    // (0..(ins_size+1)).into_iter().collect::<Vec<_>>())

    let major = major_pos_ins_vec
        .iter()
        .flat_map(|&(major_pos, ins_size)| vec![major_pos; ins_size + 1].into_iter())
        .collect::<Vec<_>>();
    let minor = major_pos_ins_vec
        .iter()
        .flat_map(|&(_, ins_size)| (0..(ins_size + 1)).into_iter())
        .collect::<Vec<_>>();

    let mut cursor = 0;
    let major_start_point = major_pos_ins_vec
        .iter()
        .map(|&(_, max_ins)| {
            let cur_point = cursor;
            cursor += max_ins + 1;
            cur_point
        })
        .collect::<Vec<_>>();
    let major_pos2major_starting_point = major_pos_ins_vec
        .iter()
        .map(|&(ma, _)| ma)
        .zip(major_start_point.into_iter())
        .collect::<HashMap<_, _>>();

    // TODO: error handling
    let first_major = major[0].clone();
    let last_major = major.last().unwrap_or(&0).clone();

    let mut msa_matrix = Array2::<u8>::from_elem((records.len(), major.len()), '-' as u8);
    records.iter().enumerate().for_each(|(idx, record)| {
        build_one_record_of_msa(
            record,
            &major_pos2major_starting_point,
            msa_matrix.slice_mut(s![idx, ..]),
        );
    });

    // println!("{:?}", msa_matrix.mapv(char::from));

    let mut plp_count = count(&msa_matrix);
    // println!("{:?}", plp_count);

    let major_depth = compute_major_depth(target_len, &records);

    plp_count
        .axis_iter_mut(Axis(1))
        .enumerate()
        .for_each(|(locus, mut plp_c)| {
            let cur_major = major[locus];
            let cur_depth = *major_depth.get(&cur_major).unwrap() as f32;
            plp_c.mapv_inplace(|v| v / cur_depth);
        });

    let mut major = major;
    let mut minor = minor;
    if first_major > 0 {
        let pad = Array2::<f32>::from_elem((4, first_major), 0.0);
        plp_count = concatenate![Axis(1), pad, plp_count];
        let mut header = (0..first_major).into_iter().collect::<Vec<usize>>();
        header.extend_from_slice(&major);
        major = header;

        let mut header = vec![0; first_major];
        header.extend_from_slice(&minor);
        minor = header;
    }
    if (last_major + 1) < target_len {
        let pad_len = target_len - last_major - 1;
        let pad = Array2::<f32>::from_elem((4, pad_len), 0.0);
        plp_count = concatenate![Axis(1), plp_count, pad];
        let tail = (last_major + 1..target_len)
            .into_iter()
            .collect::<Vec<usize>>();
        major.extend_from_slice(&tail);

        let tail = vec![0; pad_len];
        minor.extend_from_slice(&tail);
    }

    // println!("{:?}", plp_count);

    super::PlpInfo {
        normed_count: plp_count,
        major: major,
        minor: minor,
    }
}

fn build_one_record_of_msa(
    record: &Record,
    major_pos2major_starting_point: &HashMap<usize, usize>,
    mut result: ArrayViewMut1<u8>,
) {
    let record_ext = BamRecordExt::new(record);
    let ref_start = record_ext.reference_start() as i64;
    let ref_end = record_ext.reference_end() as i64;
    let q_start = record_ext.query_alignment_start() as i64;
    let q_end = record_ext.query_alignment_end() as i64;
    let query = record_ext.get_seq();
    let query = query.as_bytes();
    let mut q_pos_cursor = None;
    let mut r_pos_cursor = None;
    let mut delta = 0;
    for [qpos, rpos] in record.aligned_pairs_full() {
        if qpos.is_some() {
            q_pos_cursor = qpos;
        }
        if rpos.is_some() {
            r_pos_cursor = rpos;
        }

        if q_pos_cursor.is_none() || r_pos_cursor.is_none() {
            continue;
        }

        if q_pos_cursor.unwrap() < q_start || r_pos_cursor.unwrap() < ref_start {
            continue;
        }

        if q_pos_cursor.unwrap() >= q_end || r_pos_cursor.unwrap() >= ref_end {
            break;
        }

        if rpos.is_some() {
            delta = 0;
        } else {
            delta += 1;
        }
        let r_cursor = r_pos_cursor.map(|v| v as usize).unwrap();
        let r_cursor = &r_cursor;
        if !major_pos2major_starting_point.contains_key(r_cursor) {
            continue;
        }

        let base_pos = *major_pos2major_starting_point.get(r_cursor).unwrap();

        if let Some(qpos) = qpos.map(|v| v as usize) {
            result[base_pos + delta] = query[qpos];
        }
    }
}

pub fn compute_major_depth(target_len: usize, records: &Vec<Record>) -> HashMap<usize, usize> {
    let all_records_start_end = records
        .iter()
        .map(|record| {
            let ext = BamRecordExt::new(record);
            (ext.reference_start(), ext.reference_end())
        })
        .collect::<Vec<(usize, usize)>>();

    (0..target_len)
        .into_iter()
        .map(|major_pos| {
            let depth = all_records_start_end
                .iter()
                .filter(|(start, end)| *start <= major_pos && major_pos < *end)
                .count();
            (major_pos, depth)
        })
        .collect::<HashMap<usize, usize>>()
}

pub fn count(msa_matrix: &Array2<u8>) -> Array2<f32> {
    let mut result_matrix = Array2::<f32>::from_elem((4, msa_matrix.shape()[1]), 0.0);
    msa_matrix.axis_iter(Axis(0)).for_each(|row| {
        row.iter().enumerate().for_each(|(locus, &base)| {
            let mut base_idx = BASE2IDX[base as usize] as usize;
            if base_idx > 0 {
                base_idx -= 1;
                result_matrix[[base_idx, locus]] += 1.0;
            }
        });
    });

    result_matrix
}

#[cfg(test)]
mod test {
    use rust_htslib::bam::{Read, Record};

    use crate::pileup_counter::{extract_seq_info_from_header, plp_from_records::plp_from_records};

    #[test]
    fn test_plp_from_records() {
        let fpath = "./test-data/Group_0_Adaptor-barcode295-2.sort.bam";
        let mut reader = rust_htslib::bam::Reader::from_path(fpath).unwrap();
        reader.set_threads(40).unwrap();

        let mut record = Record::new();
        let mut records = vec![];
        loop {
            if let Some(Ok(_)) = reader.read(&mut record) {
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
            .take(10)
            .collect();
        let seq_info = extract_seq_info_from_header(reader.header()).unwrap();
        let seq_len = seq_info.length;
        plp_from_records(&records, seq_len);
    }
}
