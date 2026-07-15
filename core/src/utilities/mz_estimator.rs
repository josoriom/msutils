use std::cmp::Ordering;

use crate::utilities::{find_features::MzTolerance, math::median};

#[derive(Clone, Copy, Debug)]
pub struct SampleMz {
    pub mz: f64,
    pub intensity: f64,
}

pub trait MzEstimator {
    fn combine(&self, values: &[SampleMz], tolerance: &MzTolerance) -> Option<f64>;
}

pub struct MedianMzApex;

impl MzEstimator for MedianMzApex {
    fn combine(&self, values: &[SampleMz], tolerance: &MzTolerance) -> Option<f64> {
        let kept = keep_near_anchor(values, tolerance);
        if kept.is_empty() {
            return None;
        }
        let mut mzs: Vec<f64> = kept.iter().map(|v| v.mz).collect();
        Some(median(&mut mzs))
    }
}

pub struct WeightedIntensityMedian;

impl MzEstimator for WeightedIntensityMedian {
    fn combine(&self, values: &[SampleMz], tolerance: &MzTolerance) -> Option<f64> {
        let kept = keep_near_anchor(values, tolerance);
        weighted_median(&kept)
    }
}

pub fn highest_intensity(values: &[SampleMz]) -> Option<SampleMz> {
    values.iter().copied().max_by(|a, b| {
        a.intensity
            .partial_cmp(&b.intensity)
            .unwrap_or(Ordering::Equal)
    })
}

pub fn keep_near_anchor(values: &[SampleMz], tolerance: &MzTolerance) -> Vec<SampleMz> {
    match highest_intensity(values) {
        None => Vec::new(),
        Some(anchor) => values
            .iter()
            .copied()
            .filter(|v| tolerance.are_close_to_ref(v.mz, anchor.mz))
            .collect(),
    }
}

pub fn weighted_median(values: &[SampleMz]) -> Option<f64> {
    let total: f64 = values.iter().map(|v| v.intensity).sum();
    if !total.is_finite() || total <= 0.0 {
        return None;
    }
    let mut sorted: Vec<SampleMz> = values.to_vec();
    sorted.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap_or(Ordering::Equal));
    let half = total / 2.0;
    let mut running = 0.0;
    for value in &sorted {
        running += value.intensity;
        if running >= half {
            return Some(value.mz);
        }
    }
    sorted.last().map(|value| value.mz)
}

#[derive(Clone, Debug, Default)]
pub enum MzEstimatorKind {
    #[default]
    MedianMzApex,
    WeightedIntensityMedian,
}

pub fn make_estimator(kind: &MzEstimatorKind) -> Box<dyn MzEstimator + Send + Sync> {
    match kind {
        MzEstimatorKind::MedianMzApex => Box::new(MedianMzApex),
        MzEstimatorKind::WeightedIntensityMedian => Box::new(WeightedIntensityMedian),
    }
}

pub fn same_mass_gap(mz: f64, ppm: f64) -> f64 {
    0.5 * ppm * 1e-6 * mz
}

pub fn has_multiple_masses(centroids: &[f64], max_gap: f64) -> bool {
    let mut sorted = centroids.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    sorted.windows(2).any(|pair| pair[1] - pair[0] > max_gap)
}

pub fn split_by_gap(centroids: &[f64], max_gap: f64) -> Vec<Vec<f64>> {
    let mut sorted = centroids.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let mut groups: Vec<Vec<f64>> = Vec::new();
    for value in sorted {
        let start_new = match groups.last() {
            Some(group) => value - *group.last().unwrap() > max_gap,
            None => true,
        };
        if start_new {
            groups.push(vec![value]);
        } else {
            groups.last_mut().unwrap().push(value);
        }
    }
    groups
}
