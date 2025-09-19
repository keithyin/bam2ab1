use crate::pileup_counter::{BASE2IDX, PlpInfo};
use ndarray::{Array1, s};
use std::f32::consts::PI;

pub const MAX_PULSE: f32 = 50.;

/// peak width . for example 10
/// 1 2 3 4 5 p 6 7 8 9 10 11  12 13 14 15 16 p
/// mu means timestep
pub fn get_mu_list(length: usize, peak_width: usize) -> Vec<usize> {
    let mut mu_list = vec![];

    for idx in 0..length {
        if idx == 0 {
            mu_list.push(peak_width / 2);
        } else {
            let cur_mu = *mu_list.last().unwrap() + peak_width + 1;
            mu_list.push(cur_mu);
        }
    }

    mu_list
}

pub fn get_sigma_list(mu_list: &Vec<usize>) -> Vec<f32> {
    let k = 0.001;
    let c = 4.0;
    mu_list.iter().map(|&mu| k * (mu as f32) + c).collect()
}

pub struct Gaussian {
    mu: f32,
    sigma: f32,
    pulse: f32,

    three_sigma_law_low: f32,
    three_sigma_law_high: f32,
}

impl Gaussian {
    pub fn new(mu: f32, sigma: f32, pulse: f32) -> Self {
        let three_sigma_law_low = mu - 3.0 * sigma;
        let three_sigma_law_high = mu + 3.0 * sigma;

        Self {
            mu,
            sigma,
            pulse,
            three_sigma_law_low,
            three_sigma_law_high,
        }
    }

    pub fn compute(&self, x: f32) -> f32 {
        if x <= self.three_sigma_law_low || x >= self.three_sigma_law_high {
            return 0.0;
        }
        let v1 = self.pulse / ((2.0 * PI).sqrt() * self.sigma);
        let v2 = (-(x - self.mu).powi(2) / (2.0 * self.sigma.powi(2))).exp();
        v1 + v2
    }
}

pub fn transform_plp_info_2_ab1_data(plp_info: &PlpInfo) {
    let mu_list = get_mu_list(plp_info.major.len(), 10);
    let sigma_list = get_sigma_list(&mu_list);

    let last_mu = (*mu_list.last().unwrap()) as f32;
    let last_sigma = *sigma_list.last().unwrap();
    let last_tt = (last_mu + 3.0 * last_sigma).ceil() as usize;

    let g_pulse = plp_info
        .normed_count
        .slice(s![BASE2IDX['G' as usize] as usize, ..])
        .mapv(|v| v * 50.);

    let a_pulse = plp_info
        .normed_count
        .slice(s![BASE2IDX['A' as usize] as usize, ..])
        .mapv(|v| v * 50.);

    let t_pulse = plp_info
        .normed_count
        .slice(s![BASE2IDX['T' as usize] as usize, ..])
        .mapv(|v| v * 50.);

    let c_pulse = plp_info
        .normed_count
        .slice(s![BASE2IDX['C' as usize] as usize, ..])
        .mapv(|v| v * 50.);

    let g_data = build_data(g_pulse, &mu_list, &sigma_list, last_tt);
    let a_data = build_data(a_pulse, &mu_list, &sigma_list, last_tt);
    let t_data = build_data(t_pulse, &mu_list, &sigma_list, last_tt);
    let c_data = build_data(c_pulse, &mu_list, &sigma_list, last_tt);

    let peek_loc = mu_list
        .iter()
        .zip(plp_info.minor.iter())
        .filter(|(_, mino)| **mino == 0)
        .map(|(loc, _)| loc)
        .copied()
        .map(|v| v as u16)
        .collect::<Vec<_>>();


    

    // let
}

pub fn build_data(
    peak_pulse: Array1<f32>,
    mu_list: &Vec<usize>,
    sigma_list: &Vec<f32>,
    last_tt: usize,
) -> Vec<u16> {
    let gaussians = peak_pulse
        .iter()
        .zip(mu_list.iter())
        .zip(sigma_list.iter())
        .map(|((&pulse, &mu), &sigma)| Gaussian::new(mu as f32, sigma, pulse))
        .collect::<Vec<_>>();

    (1..last_tt)
        .into_iter()
        .map(|tt| {
            gaussians
                .iter()
                .map(|gauss| gauss.compute(tt as f32))
                .sum::<f32>() as u16
        })
        .collect()
}
