// tests/helpers.rs
use msutils::utilities::structs::{DataXY, Peak};

#[allow(dead_code)]
pub fn dump_peaks(peaks: &[Peak]) {
    println!("peaks.len() = {}", peaks.len());
    for (i, p) in peaks.iter().enumerate() {
        let r2 =
            p.r2.map(|v| format!("{:.3}", v))
                .unwrap_or_else(|| "N/A".to_string());
        println!(
            "#{:03} from={:.6} to={:.6} rt={:.6} integral={:.3} intensity={:.3} r2={} n_points={} noise={:.3}",
            i, p.from, p.to, p.rt, p.integral, p.intensity, r2, p.n_points, p.noise
        );
    }
}

#[inline]
#[allow(dead_code)]
pub fn gaussian_value(x: f64, mu: f64, sigma: f64, amp: f64, base: f64) -> f64 {
    base + amp * (-0.5 * ((x - mu) / sigma).powi(2)).exp()
}

#[allow(dead_code)]
pub fn gaussian_mixture_f32(
    xs: &[f64],
    peaks: &[(f64, f64, f64)],
    base: f64,
    noise: f64,
) -> Vec<f32> {
    xs.iter()
        .map(|&x| {
            let mut y = base;
            for &(mu, sigma, amp) in peaks {
                y += gaussian_value(x, mu, sigma, amp, 0.0);
            }
            if noise > 0.0 {
                let z = ((x * 137.13).sin() + (x * 73.7).cos()) * 0.5;
                y += z * noise;
            }
            y as f32
        })
        .collect()
}

pub fn make_grid(start: f64, end: f64, n: usize) -> Vec<f64> {
    if n <= 1 {
        return vec![start];
    }
    (0..n)
        .map(|i| start + (end - start) * (i as f64) / ((n - 1) as f64))
        .collect()
}

#[allow(dead_code)]
pub fn linspace(from: f64, to: f64, n: usize) -> Vec<f64> {
    make_grid(from, to, n)
}

#[allow(dead_code)]
pub fn jitter(i: u32) -> f64 {
    let mut x = i.wrapping_mul(1664525).wrapping_add(1013904223);
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    (x as f64 / (u32::MAX as f64)) - 0.5
}

#[allow(dead_code)]
pub fn approx_eq(a: f64, b: f64, tol: f64) -> bool {
    (a - b).abs() <= tol
}

#[allow(dead_code)]
pub fn data_xy(xs: Vec<f64>, ys: Vec<f64>) -> DataXY {
    DataXY { x: xs, y: ys }
}

#[allow(dead_code)]
pub fn uniform_vec_f32(n: usize, lo: f32, hi: f32, seed: u64) -> Vec<f32> {
    assert!(hi > lo);
    let mut out = Vec::with_capacity(n);
    let mut s = seed | 1; // odd
    for _ in 0..n {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let u = ((s >> 11) as f64) * (1.0 / (1u64 << 53) as f64);
        let v = lo as f64 + (hi - lo) as f64 * u;
        out.push(v as f32);
    }
    out
}

#[allow(dead_code)]
pub fn shuffle_with_seed<T>(xs: &mut [T], seed: u64) {
    let mut s = seed | 1;
    let n = xs.len();
    if n <= 1 {
        return;
    }
    for i in (1..n).rev() {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        let u = ((s >> 11) as f64) * (1.0 / (1u64 << 53) as f64);
        let j = (u * ((i + 1) as f64)).floor() as usize;
        xs.swap(i, j);
    }
}

use std::{collections::BTreeMap, path::Path};

use ionic::{
    ion::{IonReader, ReadOptions},
    mzml::structs::{Chromatogram, NumericArray},
};

#[allow(dead_code)]
pub struct Eic {
    pub time: Vec<f64>,
    pub intensity: Vec<f64>,
    pub params: BTreeMap<String, String>,
}

#[allow(dead_code)]
pub fn load_chromatograms(file_name: &str) -> BTreeMap<String, Eic> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join(file_name);
    load_chromatograms_from(&path)
}

#[allow(dead_code)]
pub fn load_chromatograms_from(path: &Path) -> BTreeMap<String, Eic> {
    let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("cannot read {:?}: {}", path, e));
    let mut reader = IonReader::open(
        &bytes,
        ReadOptions {
            verify_checksums: false, // TODO: update the fixtures to the lastest version of Ionic
            ..Default::default()
        },
    )
    .expect("open ion fixture");
    let mzml = reader.to_mzml().expect("decode ion fixture");

    let mut out = BTreeMap::new();
    if let Some(list) = mzml.run.chromatogram_list {
        for item in list.chromatograms {
            let time = chromatogram_array(&item, "MS:1000595");
            let intensity = chromatogram_array(&item, "MS:1000515");
            let params = item
                .user_params
                .into_iter()
                .filter_map(|param| param.value.map(|value| (param.name, value)))
                .collect();
            out.insert(
                item.id,
                Eic {
                    time,
                    intensity,
                    params,
                },
            );
        }
    }
    out
}

#[allow(dead_code)]
fn chromatogram_array(item: &Chromatogram, accession: &str) -> Vec<f64> {
    let arrays = &item
        .binary_data_array_list
        .as_ref()
        .expect("binary arrays")
        .binary_data_arrays;
    for array in arrays {
        if array
            .cv_params
            .iter()
            .any(|p| p.accession.as_deref() == Some(accession))
        {
            return numeric_to_f64(array.binary.as_ref().expect("binary payload"));
        }
    }
    panic!("array {accession} missing in {}", item.id);
}

#[allow(dead_code)]
fn numeric_to_f64(array: &NumericArray) -> Vec<f64> {
    match array {
        NumericArray::F64(v) => v.clone(),
        NumericArray::F32(v) => v.iter().map(|&x| x as f64).collect(),
        NumericArray::F16(v) => v.iter().map(|&x| x as f64).collect(),
        NumericArray::I64(v) => v.iter().map(|&x| x as f64).collect(),
        NumericArray::I32(v) => v.iter().map(|&x| x as f64).collect(),
        NumericArray::I16(v) => v.iter().map(|&x| x as f64).collect(),
    }
}

#[allow(dead_code)]
pub fn group_ids<'a>(chromatograms: &'a BTreeMap<String, Eic>, prefix: &str) -> Vec<&'a String> {
    chromatograms
        .keys()
        .filter(|id| id.starts_with(prefix))
        .collect()
}

#[allow(dead_code)]
pub fn param_f64(eic: &Eic, key: &str) -> f64 {
    eic.params
        .get(key)
        .unwrap_or_else(|| panic!("param {key} missing"))
        .parse()
        .expect("numeric param")
}

#[allow(dead_code)]
pub fn param_bool(eic: &Eic, key: &str) -> bool {
    eic.params
        .get(key)
        .map(|value| value == "true")
        .unwrap_or(false)
}
