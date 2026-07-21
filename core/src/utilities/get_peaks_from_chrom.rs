use ionic::mzml::structs::{
    BinaryDataArray, BinaryDataArrayList, Chromatogram, MzML, NumericArray,
};
#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use rayon::{ThreadPoolBuilder, prelude::*};

use crate::utilities::{
    calculate_baseline::BaselineOptions,
    find_peaks::{FindPeaksOptions, find_peaks},
    get_peak::get_peak,
    structs::{DataXY, Roi},
};

const ACC_TIME_ARRAY: &str = "MS:1000595";
const ACC_INTENSITY_ARRAY: &str = "MS:1000515";

#[derive(Clone, Debug)]
pub struct ChromPeakRow {
    pub index: usize,
    pub id: String,
    pub target_rt: f64,
    pub peak_rt: f64,
    pub from_rt: f64,
    pub to_rt: f64,
    pub intensity: f64,
    pub area: f64,
    pub total_area: f64,
    pub timestamp: String,
}

impl ChromPeakRow {
    fn empty(index: usize, id: String, target_rt: f64, timestamp: String) -> Self {
        Self {
            index,
            id,
            target_rt,
            peak_rt: 0.0,
            from_rt: 0.0,
            to_rt: 0.0,
            intensity: 0.0,
            area: 0.0,
            total_area: 0.0,
            timestamp,
        }
    }
}

pub fn get_peaks_from_chrom(
    mzml: &MzML,
    items: &[Roi],
    options: Option<FindPeaksOptions>,
    cores: usize,
) -> Option<Vec<ChromPeakRow>> {
    let run = &mzml.run;
    let chroms = run.chromatogram_list.as_ref()?.chromatograms.as_slice();
    let timestamp = run.start_time_stamp.clone().unwrap_or_default();

    let mut opts = options.unwrap_or_default();
    opts.baseline = Some(BaselineOptions {
        lambda: Some(50.0),
        ..Default::default()
    });
    let opts = Some(opts);

    let make_row = |roi: &Roi| {
        let timestamp = timestamp.clone();
        let Roi::Chrom {
            id,
            sample_index,
            rt,
            range,
        } = roi
        else {
            return ChromPeakRow::empty(0, "invalid".to_string(), 0.0, timestamp);
        };
        let (item_index, rt, range) = (*sample_index, *rt, *range);

        if range <= 0.0 || !rt.is_finite() {
            return ChromPeakRow::empty(item_index, id.clone(), rt, timestamp);
        }

        if item_index >= chroms.len() {
            return ChromPeakRow::empty(item_index, id.clone(), rt, timestamp);
        }

        let chrom = &chroms[item_index];

        let Some((time_values, intensity_values)) = chromatogram_xy(chrom) else {
            let chrom_index = chrom.index.unwrap_or(item_index as u32) as usize;
            return ChromPeakRow::empty(chrom_index, chrom.id.clone(), rt, timestamp);
        };

        let chrom_index = chrom.index.unwrap_or(item_index as u32) as usize;
        compute_one(
            chrom_index,
            chrom.id.as_str(),
            time_values,
            intensity_values,
            roi,
            timestamp,
            &opts,
        )
    };

    #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
    {
        return Some(items.iter().map(make_row).collect());
    }

    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    {
        if cores <= 1 || items.len() < 2 {
            return Some(items.iter().map(make_row).collect());
        }
        let pool = ThreadPoolBuilder::new().num_threads(cores).build().ok()?;
        Some(pool.install(|| items.par_iter().map(make_row).collect()))
    }
}

fn chromatogram_xy(ch: &Chromatogram) -> Option<(Vec<f64>, Vec<f32>)> {
    let bal = ch.binary_data_array_list.as_ref()?;

    let time_bda =
        find_bda_by_accession(bal, ACC_TIME_ARRAY).or_else(|| bal.binary_data_arrays.first());
    let inten_bda =
        find_bda_by_accession(bal, ACC_INTENSITY_ARRAY).or_else(|| bal.binary_data_arrays.get(1));

    let (time_bda, inten_bda) = (time_bda?, inten_bda?);

    let mut x = bda_to_f64(time_bda);
    let mut y = bda_to_f32(inten_bda);

    let n = x.len().min(y.len());
    x.truncate(n);
    y.truncate(n);

    if n >= 1 { Some((x, y)) } else { None }
}

#[inline]
fn find_bda_by_accession<'a>(
    bal: &'a BinaryDataArrayList,
    accession: &str,
) -> Option<&'a BinaryDataArray> {
    bal.binary_data_arrays.iter().find(|bda| {
        bda.cv_params
            .iter()
            .any(|cv| cv.accession.as_deref() == Some(accession))
    })
}

#[inline]
fn bda_to_f64(bda: &BinaryDataArray) -> Vec<f64> {
    let Some(bin) = bda.binary.as_ref() else {
        return Vec::new();
    };
    match bin {
        NumericArray::F64(v) => v.clone(),
        NumericArray::F32(v) => v.iter().map(|&x| x as f64).collect(),
        NumericArray::F16(v) => v.iter().map(|&x| x as f64).collect(),
        NumericArray::I64(v) => v.iter().map(|&x| x as f64).collect(),
        NumericArray::I32(v) => v.iter().map(|&x| x as f64).collect(),
        NumericArray::I16(v) => v.iter().map(|&x| x as f64).collect(),
    }
}

#[inline]
fn bda_to_f32(bda: &BinaryDataArray) -> Vec<f32> {
    let Some(bin) = bda.binary.as_ref() else {
        return Vec::new();
    };
    match bin {
        NumericArray::F16(v) => v.iter().map(|&x| x as f32).collect(),
        NumericArray::F32(v) => v.clone(),
        NumericArray::F64(v) => v.iter().map(|&x| x as f32).collect(),
        NumericArray::I64(v) => v.iter().map(|&x| x as f32).collect(),
        NumericArray::I32(v) => v.iter().map(|&x| x as f32).collect(),
        NumericArray::I16(v) => v.iter().map(|&x| x as f32).collect(),
    }
}

#[inline]
fn compute_one(
    chrom_index: usize,
    chrom_id: &str,
    time_values: Vec<f64>,
    intensity_values: Vec<f32>,
    roi: &Roi,
    timestamp: String,
    options: &Option<FindPeaksOptions>,
) -> ChromPeakRow {
    let mut peak_rt = 0.0_f64;
    let mut from_rt = 0.0_f64;
    let mut to_rt = 0.0_f64;
    let mut intensity = 0.0_f64;
    let mut area = 0.0_f64;
    let mut total_area = 0.0_f64;
    let Roi::Chrom { rt, range, .. } = roi else {
        return ChromPeakRow::empty(0, "invalid".to_string(), 0.0, timestamp);
    };

    if time_values.len() >= 3 && time_values.len() == intensity_values.len() {
        let mut signal_values: Vec<f64> = Vec::with_capacity(intensity_values.len());
        for &value in &intensity_values {
            signal_values.push(value as f64);
        }

        let data = DataXY {
            x: time_values,
            y: signal_values,
        };

        let picked = get_peak(&data, &Roi::peak(*rt, *range), options.clone());
        if picked.intensity > 0.0 {
            peak_rt = picked.rt;
            from_rt = picked.from;
            to_rt = picked.to;
            intensity = picked.intensity;
            area = picked.integral;
        }

        let peaks = find_peaks(&data, options.clone());
        for peak in peaks {
            total_area += peak.integral;
        }
    }

    ChromPeakRow {
        index: chrom_index,
        id: chrom_id.trim_end().to_string(),
        target_rt: *rt,
        peak_rt,
        from_rt,
        to_rt,
        intensity,
        area,
        total_area,
        timestamp,
    }
}
