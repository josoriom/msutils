use std::collections::HashMap;

use ionic::ScanSource;
use serde::Serialize;

use super::calculate_eic::summed_intensity_in_window;

#[derive(Clone, Debug, Serialize)]
pub struct IonImage {
    pub width: u32,
    pub height: u32,
    pub min_x: u32,
    pub min_y: u32,
    pub min_z: u32,
    pub max_z: u32,
    pub data: Vec<f64>,
    pub counts: Vec<u32>,
}

impl IonImage {
    fn empty() -> Self {
        Self {
            width: 0,
            height: 0,
            min_x: 0,
            min_y: 0,
            min_z: 0,
            max_z: 0,
            data: Vec::new(),
            counts: Vec::new(),
        }
    }
}

pub fn compute_ion_image(
    source: &mut impl ScanSource,
    target_mz: f64,
    tolerance: f64,
    ms_level: u8,
) -> IonImage {
    let lower = target_mz - tolerance;
    let upper = target_mz + tolerance;

    let mut candidates: Vec<(usize, u32, u32, u32)> = Vec::new();
    source.for_each_summary(&mut |index, summary| {
        if ms_level != 0 && summary.ms_level != ms_level {
            return;
        }
        candidates.push((
            index,
            summary.position_x,
            summary.position_y,
            summary.position_z,
        ));
    });

    let mut accumulator: HashMap<(u32, u32), (f64, u32)> = HashMap::new();
    let mut min_x = u32::MAX;
    let mut max_x = 0u32;
    let mut min_y = u32::MAX;
    let mut max_y = 0u32;
    let mut min_z = u32::MAX;
    let mut max_z = 0u32;

    let mut mz_buffer = Vec::new();
    let mut intensity_buffer = Vec::new();
    for (index, x, y, z) in &candidates {
        if !source.load_scan(*index, &mut mz_buffer, &mut intensity_buffer) {
            continue;
        }
        let len = mz_buffer.len().min(intensity_buffer.len());
        let sum =
            summed_intensity_in_window(&mz_buffer[..len], &intensity_buffer[..len], lower, upper);
        let entry = accumulator.entry((*x, *y)).or_insert((0.0, 0));
        entry.0 += sum;
        entry.1 += 1;
        if *x < min_x {
            min_x = *x;
        }
        if *x > max_x {
            max_x = *x;
        }
        if *y < min_y {
            min_y = *y;
        }
        if *y > max_y {
            max_y = *y;
        }
        if *z < min_z {
            min_z = *z;
        }
        if *z > max_z {
            max_z = *z;
        }
    }

    if accumulator.is_empty() {
        return IonImage::empty();
    }

    let width = max_x - min_x + 1;
    let height = max_y - min_y + 1;
    let cells = (width as usize) * (height as usize);
    let mut data = vec![0.0f64; cells];
    let mut counts = vec![0u32; cells];
    for ((x, y), (sum, count)) in accumulator {
        let cell = ((y - min_y) as usize) * (width as usize) + (x - min_x) as usize;
        counts[cell] = count;
        data[cell] = if count > 0 { sum / count as f64 } else { 0.0 };
    }

    IonImage {
        width,
        height,
        min_x,
        min_y,
        min_z,
        max_z,
        data,
        counts,
    }
}
