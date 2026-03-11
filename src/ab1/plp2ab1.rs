use crate::{ab1::ab1_file, pileup_counter::PlpInfo};

pub trait TPlp2Ab1 {
    fn transform(
        &self,
        plp_info: &PlpInfo,
        target_seq: &str,
        peak_width: Option<usize>,
        sample_name: Option<String>,
    ) -> ab1_file::AbiFile;
}
