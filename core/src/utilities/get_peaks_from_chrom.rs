use rayon::{ThreadPoolBuilder, prelude::*};

use crate::utilities::{
    calculate_baseline::BaselineOptions,
    find_peaks::{FindPeaksOptions, find_peaks},
    get_peak::get_peak,
    structs::{ChromRoi, DataXY, Roi},
};

use octo::{BinaryData, BinaryDataArray, BinaryDataArrayList, Chromatogram, MzML};

const ACC_TIME_ARRAY: &str = "MS:1000595";
const ACC_INTENSITY_ARRAY: &str = "MS:1000515";

pub fn get_peaks_from_chrom(
    mzml: &MzML,
    items: &[ChromRoi],
    options: Option<FindPeaksOptions>,
    cores: usize,
) -> Option<Vec<(usize, String, f64, f64, f64, f64, f64, f64, f64, String)>> {
    let run = &mzml.run;
    let chroms = run.chromatogram_list.as_ref()?.chromatograms.as_slice();
    let time_stamp = run.start_time_stamp.clone().unwrap_or_default();

    let mut opts = options.unwrap_or_default();
    opts.baseline_options = Some(BaselineOptions {
        baseline_window: Some(50.0),
        ..Default::default()
    });
    let opts = Some(opts);

    let f = |roi: &ChromRoi| {
        let ts = time_stamp.clone();

        if roi.window <= 0.0 || !roi.rt.is_finite() {
            return (
                roi.idx,
                roi.id.clone(),
                roi.rt,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                ts,
            );
        }

        let i = roi.idx;
        if i >= chroms.len() {
            return (i, roi.id.clone(), roi.rt, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ts);
        }

        let ch = &chroms[i];

        let Some((x, y)) = chromatogram_xy(ch) else {
            let ch_index = ch.index.unwrap_or(i as u32) as usize;
            return (
                ch_index,
                ch.id.clone(),
                roi.rt,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                ts,
            );
        };

        let ch_index = ch.index.unwrap_or(i as u32) as usize;
        compute_one(ch_index, ch.id.as_str(), x, y, roi, ts, &opts)
    };

    if cores <= 1 || items.len() < 2 {
        Some(items.iter().map(f).collect())
    } else {
        let pool = ThreadPoolBuilder::new().num_threads(cores).build().ok()?;
        Some(pool.install(|| items.par_iter().map(f).collect()))
    }
}

fn chromatogram_xy(ch: &Chromatogram) -> Option<(Vec<f64>, Vec<f32>)> {
    let bal = ch.binary_data_array_list.as_ref()?;

    let time_bda =
        find_bda_by_accession(bal, ACC_TIME_ARRAY).or_else(|| bal.binary_data_arrays.get(0));
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
        BinaryData::F64(v) => v.clone(),
        BinaryData::F32(v) => v.iter().map(|&x| x as f64).collect(),
        BinaryData::F16(v) => v.iter().map(|&x| x as f64).collect(),
        BinaryData::I64(v) => v.iter().map(|&x| x as f64).collect(),
        BinaryData::I32(v) => v.iter().map(|&x| x as f64).collect(),
        BinaryData::I16(v) => v.iter().map(|&x| x as f64).collect(),
    }
}

#[inline]
fn bda_to_f32(bda: &BinaryDataArray) -> Vec<f32> {
    let Some(bin) = bda.binary.as_ref() else {
        return Vec::new();
    };

    match bin {
        BinaryData::F16(v) => v.iter().map(|&x| x as f32).collect(),
        BinaryData::F32(v) => v.clone(),
        BinaryData::F64(v) => v.iter().map(|&x| x as f32).collect(),
        BinaryData::I64(v) => v.iter().map(|&x| x as f32).collect(),
        BinaryData::I32(v) => v.iter().map(|&x| x as f32).collect(),
        BinaryData::I16(v) => v.iter().map(|&x| x as f32).collect(),
    }
}

#[inline]
fn compute_one(
    ch_index: usize,
    ch_id: &str,
    x: Vec<f64>,
    y: Vec<f32>,
    roi: &ChromRoi,
    time_stamp: String,
    options: &Option<FindPeaksOptions>,
) -> (usize, String, f64, f64, f64, f64, f64, f64, f64, String) {
    let mut rt = 0.0_f64;
    let mut from = 0.0_f64;
    let mut to = 0.0_f64;
    let mut intensity_val = 0.0_f64;
    let mut integral = 0.0_f64;
    let mut total_area = 0.0_f64;

    if x.len() >= 3 && x.len() == y.len() {
        let mut y64: Vec<f64> = Vec::with_capacity(y.len());
        for &v in &y {
            y64.push(v as f64);
        }

        let data = DataXY { x, y: y64 };

        if let Some(p) = get_peak(
            &data,
            Roi {
                rt: roi.rt,
                window: roi.window,
            },
            options.clone(),
        ) {
            rt = p.rt;
            from = p.from;
            to = p.to;
            intensity_val = p.intensity;
            integral = p.integral;
        }

        let peaks = find_peaks(&data, options.clone());
        for peak in peaks {
            total_area += peak.integral;
        }
    }

    (
        ch_index,
        ch_id.trim_end().to_string(),
        roi.rt,
        rt,
        from,
        to,
        intensity_val,
        integral,
        total_area,
        time_stamp,
    )
}
