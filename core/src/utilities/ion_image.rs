use std::collections::HashMap;

use ionic::ion::{ByteRange, IonReader, Range};
use ionic::{ScanSource, ScanSummary};
use serde::Serialize;

use super::calculate_eic::{
    EicReader, FastError, read_mz_window, sort_and_dedup_ranges, summed_intensity_in_window,
};

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

type Candidate = (usize, u32, u32, u32);
type PixelSums = HashMap<(u32, u32), (f64, u32)>;

struct PixelBounds {
    seen: bool,
    min_x: u32,
    max_x: u32,
    min_y: u32,
    max_y: u32,
    min_z: u32,
    max_z: u32,
}

impl PixelBounds {
    fn new() -> Self {
        Self {
            seen: false,
            min_x: u32::MAX,
            max_x: 0,
            min_y: u32::MAX,
            max_y: 0,
            min_z: u32::MAX,
            max_z: 0,
        }
    }

    fn update(&mut self, x: u32, y: u32, z: u32) {
        self.seen = true;
        if x < self.min_x {
            self.min_x = x;
        }
        if x > self.max_x {
            self.max_x = x;
        }
        if y < self.min_y {
            self.min_y = y;
        }
        if y > self.max_y {
            self.max_y = y;
        }
        if z < self.min_z {
            self.min_z = z;
        }
        if z > self.max_z {
            self.max_z = z;
        }
    }
}

fn scan_is_wanted(scan_level: u8, wanted_level: u8) -> bool {
    wanted_level == 0 || scan_level == wanted_level
}

fn collect_candidates(reader: &mut EicReader, ms_level: u8) -> Vec<Candidate> {
    let mut candidates: Vec<Candidate> = Vec::new();
    let mut collect = |index: usize, summary: ScanSummary| {
        if scan_is_wanted(summary.ms_level, ms_level) {
            candidates.push((
                index,
                summary.position_x,
                summary.position_y,
                summary.position_z,
            ));
        }
    };
    match reader {
        EicReader::Ion(ion) => ion.for_each_summary(&mut collect),
        EicReader::Mzml(mzml) => mzml.for_each_summary(&mut collect),
    }
    candidates
}

fn ranges_for_candidates(
    ion: &mut IonReader,
    candidates: &[Candidate],
    lower: f64,
    upper: f64,
) -> Result<Vec<ByteRange>, FastError> {
    ion.require_bounds().map_err(FastError::from)?;
    let mut ranges = Vec::new();
    for candidate in candidates {
        let scan_ranges = ion
            .byte_ranges(candidate.0, Range { from: lower, to: upper })
            .map_err(FastError::from)?;
        ranges.extend(scan_ranges);
    }
    sort_and_dedup_ranges(&mut ranges);
    Ok(ranges)
}

fn fold_candidates(
    reader: &mut EicReader,
    candidates: &[Candidate],
    lower: f64,
    upper: f64,
    sums: &mut PixelSums,
    bounds: &mut PixelBounds,
) {
    let mut mz_buffer = Vec::new();
    let mut intensity_buffer = Vec::new();
    for (index, x, y, z) in candidates {
        if read_mz_window(reader, *index, lower, upper, &mut mz_buffer, &mut intensity_buffer)
            .is_err()
        {
            continue;
        }
        let len = mz_buffer.len().min(intensity_buffer.len());
        let sum =
            summed_intensity_in_window(&mz_buffer[..len], &intensity_buffer[..len], lower, upper);
        let entry = sums.entry((*x, *y)).or_insert((0.0, 0));
        entry.0 += sum;
        entry.1 += 1;
        bounds.update(*x, *y, *z);
    }
}

fn build_image(sums: &PixelSums, bounds: &PixelBounds) -> IonImage {
    if !bounds.seen || sums.is_empty() {
        return IonImage::empty();
    }

    let width = bounds.max_x - bounds.min_x + 1;
    let height = bounds.max_y - bounds.min_y + 1;
    let cells = (width as usize) * (height as usize);
    let mut data = vec![0.0f64; cells];
    let mut counts = vec![0u32; cells];
    for ((x, y), (sum, count)) in sums {
        let cell = ((y - bounds.min_y) as usize) * (width as usize) + (x - bounds.min_x) as usize;
        counts[cell] = *count;
        data[cell] = if *count > 0 { sum / *count as f64 } else { 0.0 };
    }

    IonImage {
        width,
        height,
        min_x: bounds.min_x,
        min_y: bounds.min_y,
        min_z: bounds.min_z,
        max_z: bounds.max_z,
        data,
        counts,
    }
}

pub fn plan_ion_image_ranges(
    ion: &mut IonReader,
    target_mz: f64,
    tolerance: f64,
    ms_level: u8,
) -> Result<Vec<ByteRange>, FastError> {
    if !target_mz.is_finite() || target_mz <= 0.0 {
        return Err(FastError::InvalidRequest);
    }
    if !tolerance.is_finite() || tolerance < 0.0 {
        return Err(FastError::InvalidRequest);
    }

    let lower = target_mz - tolerance;
    let upper = target_mz + tolerance;

    let mut candidates: Vec<Candidate> = Vec::new();
    ion.for_each_summary(&mut |index, summary| {
        if scan_is_wanted(summary.ms_level, ms_level) {
            candidates.push((index, 0, 0, 0));
        }
    });

    ranges_for_candidates(ion, &candidates, lower, upper)
}

pub fn compute_ion_image(
    reader: &mut EicReader,
    target_mz: f64,
    tolerance: f64,
    ms_level: u8,
) -> IonImage {
    let lower = target_mz - tolerance;
    let upper = target_mz + tolerance;

    let candidates = collect_candidates(reader, ms_level);
    let mut sums = PixelSums::new();
    let mut bounds = PixelBounds::new();
    fold_candidates(reader, &candidates, lower, upper, &mut sums, &mut bounds);
    build_image(&sums, &bounds)
}

pub struct ImageSession {
    target_mz: f64,
    tolerance: f64,
    candidates: Vec<Candidate>,
    sums: PixelSums,
    bounds: PixelBounds,
}

pub fn image_session_begin(
    reader: &mut EicReader,
    target_mz: f64,
    tolerance: f64,
    ms_level: u8,
) -> ImageSession {
    ImageSession {
        target_mz,
        tolerance,
        candidates: collect_candidates(reader, ms_level),
        sums: PixelSums::new(),
        bounds: PixelBounds::new(),
    }
}

impl ImageSession {
    pub fn scan_count(&self) -> usize {
        self.candidates.len()
    }

    fn batch(&self, from: usize, count: usize) -> &[Candidate] {
        let start = from.min(self.candidates.len());
        let end = (from + count).min(self.candidates.len());
        &self.candidates[start..end]
    }

    pub fn ranges(
        &self,
        ion: &mut IonReader,
        from: usize,
        count: usize,
    ) -> Result<Vec<ByteRange>, FastError> {
        let lower = self.target_mz - self.tolerance;
        let upper = self.target_mz + self.tolerance;
        ranges_for_candidates(ion, self.batch(from, count), lower, upper)
    }

    pub fn fold(&mut self, reader: &mut EicReader, from: usize, count: usize) {
        let lower = self.target_mz - self.tolerance;
        let upper = self.target_mz + self.tolerance;
        let start = from.min(self.candidates.len());
        let end = (from + count).min(self.candidates.len());
        let batch = &self.candidates[start..end];
        fold_candidates(reader, batch, lower, upper, &mut self.sums, &mut self.bounds);
    }

    pub fn finish(&self) -> IonImage {
        build_image(&self.sums, &self.bounds)
    }
}
