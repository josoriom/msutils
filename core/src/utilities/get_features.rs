use std::{
    cmp::Ordering,
    fmt::{Display, Formatter},
    io::Error,
    sync::Arc,
};

use serde::Serialize;

use crate::utilities::structs::ser_finite_f64;

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use crate::utilities::gpu::GpuContext;

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use crate::utilities::parallel::run_with_cores;

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use rayon::prelude::*;

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use ionic::ion::{IonReader, ReadOptions};

use crate::utilities::{
    calculate_eic::{
        CentroidScan, EicOptions, EicReader, MS1_LEVEL, get_scan_times, lower_bound,
        mz_tolerance_for, read_mz_window, upper_bound,
    },
    find_features::{Feature, FeatureError, FindFeaturesOptions, MzTolerance, find_features},
    find_peaks::FindPeaksOptions,
    get_peak::get_peak,
    structs::{DataXY, FromTo, Roi},
};

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
const ION_CACHE_BYTES: usize = 128 * 1024 * 1024;

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use ionic::{mzml::structs::MzML, parse_mzml};
#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use memmap2::Mmap;
#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use std::{
    fs,
    fs::File,
    path::{Path, PathBuf},
    time::Instant,
};

#[derive(Debug)]
pub enum AlignmentError {
    Io(Error),
    Parse {
        path: String,
        source: String,
    },
    UnsupportedFormat(String),
    NoFiles,
    FeatureDetection {
        path: String,
        source: FeatureError,
    },
    FastPath {
        path: String,
        source: crate::utilities::calculate_eic::FastError,
    },
}

impl Display for AlignmentError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(e) => write!(f, "IO error: {}", e),
            Self::Parse { path, source } => write!(f, "parse error for {}: {}", path, source),
            Self::UnsupportedFormat(s) => write!(f, "unsupported file format: {}", s),
            Self::NoFiles => write!(f, "no valid files found in directory"),
            Self::FeatureDetection { path, source } => {
                write!(f, "feature detection error in {}: {}", path, source)
            }
            Self::FastPath { path, source } => {
                write!(f, "fast path error in {}: {}", path, source)
            }
        }
    }
}

impl std::error::Error for AlignmentError {}

impl From<Error> for AlignmentError {
    fn from(e: Error) -> Self {
        Self::Io(e)
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct MissingSample {
    pub sample: String,
    #[serde(serialize_with = "ser_finite_f64")]
    pub target_mz: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    pub center_rt: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    pub rt_from: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    pub rt_to: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct ConsensusFeature {
    #[serde(serialize_with = "ser_finite_f64")]
    pub mz: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    pub rt: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    pub from: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    pub to: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    pub intensity: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    pub integral: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    pub frequency: f64,
    #[serde(skip)]
    pub n_samples: usize,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub debug_missing: Vec<MissingSample>,
}

#[derive(Clone, Debug)]
pub struct AlignmentOptions {
    pub mz_tolerance: MzTolerance,
    pub rt_tolerance: f64,
    pub min_samples: usize,
    pub eic_options: EicOptions,
    pub peak_options: Option<FindPeaksOptions>,
}

impl Default for AlignmentOptions {
    fn default() -> Self {
        Self {
            mz_tolerance: MzTolerance {
                mz_absolute: 0.005,
                ppm: 20.0,
            },
            rt_tolerance: 0.05,
            min_samples: 1,
            eic_options: EicOptions {
                ppm_tolerance: 20.0,
                mz_tolerance: 0.005,
                ..Default::default()
            },
            peak_options: None,
        }
    }
}

#[derive(Debug)]
pub struct TaggedFeature {
    pub sample_index: usize,
    pub feature: Feature,
}

pub struct FeatureClusterer {
    pub tolerance: MzTolerance,
    pub rt_tolerance: f64,
}

type Cluster = Vec<TaggedFeature>;

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
pub enum SampleSourceKind {
    Mzml(PathBuf),
    Ion(PathBuf),
}
#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
type SampleDataset = (String, SampleSourceKind, Vec<Feature>);

pub(crate) struct SearchBounds {
    pub(crate) target_mz: f64,
    pub(crate) center_rt: f64,
    pub(crate) rt_from: f64,
    pub(crate) rt_to: f64,
}

#[derive(Debug)]
pub(crate) struct MzRtCluster {
    items: Vec<TaggedFeature>,
    sorted_mzs: Vec<f64>,
    sorted_rts: Vec<f64>,
    pub(crate) cached_median_mz: f64,
    pub(crate) cached_median_rt: f64,
}

impl MzRtCluster {
    pub(crate) fn new(item: TaggedFeature) -> Self {
        let mz = item.feature.mz;
        let rt = item.feature.rt;
        Self {
            sorted_mzs: vec![mz],
            sorted_rts: vec![rt],
            items: vec![item],
            cached_median_mz: mz,
            cached_median_rt: rt,
        }
    }

    pub(crate) fn push(&mut self, item: TaggedFeature) {
        let mz = item.feature.mz;
        let rt = item.feature.rt;
        let mz_pos = self.sorted_mzs.partition_point(|&v| v < mz);
        let rt_pos = self.sorted_rts.partition_point(|&v| v < rt);
        self.sorted_mzs.insert(mz_pos, mz);
        self.sorted_rts.insert(rt_pos, rt);
        self.items.push(item);
        self.cached_median_mz = self.sorted_mzs[self.sorted_mzs.len() / 2];
        self.cached_median_rt = self.sorted_rts[self.sorted_rts.len() / 2];
    }

    fn into_items(self) -> Vec<TaggedFeature> {
        self.items
    }
}

impl FeatureClusterer {
    pub(crate) fn cluster(&self, mut tagged: Vec<TaggedFeature>) -> Vec<Cluster> {
        tagged.sort_unstable_by(|a, b| {
            a.feature
                .mz
                .partial_cmp(&b.feature.mz)
                .unwrap_or(Ordering::Equal)
        });
        let mz_groups = self.group_by_mz(tagged);
        mz_groups
            .into_iter()
            .flat_map(|group| self.subdivide_by_rt(group.into_items()))
            .filter_map(|growing| {
                let med_mz = growing.cached_median_mz;
                let med_rt = growing.cached_median_rt;
                let kept: Vec<TaggedFeature> = growing
                    .into_items()
                    .into_iter()
                    .filter(|t| {
                        self.tolerance.are_close_to_ref(t.feature.mz, med_mz)
                            && (t.feature.rt - med_rt).abs() <= self.rt_tolerance
                    })
                    .collect();
                if kept.is_empty() { None } else { Some(kept) }
            })
            .collect()
    }

    fn group_by_mz(&self, items: Vec<TaggedFeature>) -> Vec<MzRtCluster> {
        let mut groups: Vec<MzRtCluster> = Vec::new();
        for item in items {
            let belongs = groups.last().is_some_and(|g| {
                self.tolerance
                    .are_close_to_ref(item.feature.mz, g.cached_median_mz)
            });
            if belongs {
                groups.last_mut().unwrap().push(item);
            } else {
                groups.push(MzRtCluster::new(item));
            }
        }
        groups
    }

    fn subdivide_by_rt(&self, group: Vec<TaggedFeature>) -> Vec<MzRtCluster> {
        let mut rt_sorted = group;
        rt_sorted.sort_unstable_by(|a, b| {
            a.feature
                .rt
                .partial_cmp(&b.feature.rt)
                .unwrap_or(Ordering::Equal)
        });
        let mut clusters: Vec<MzRtCluster> = Vec::new();
        for item in rt_sorted {
            let belongs = clusters
                .last()
                .is_some_and(|g| (item.feature.rt - g.cached_median_rt).abs() <= self.rt_tolerance);
            if belongs {
                clusters.last_mut().unwrap().push(item);
            } else {
                clusters.push(MzRtCluster::new(item));
            }
        }
        clusters
    }
}

type ClusterSlot = (Vec<Option<Feature>>, SearchBounds);

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
pub fn get_features(
    directory_path: &str,
    time_window: FromTo,
    feature_config: FindFeaturesOptions,
    alignment_config: AlignmentOptions,
    cores: usize,
) -> Result<Vec<ConsensusFeature>, AlignmentError> {
    let datasets = load_sample_files(directory_path)?;

    let mut feature_config = feature_config;
    if feature_config.use_gpu && feature_config.gpu_context.is_none() {
        feature_config.gpu_context = GpuContext::try_init().map(Arc::new);
    }

    let mut datasets = detect_features_per_sample(datasets, time_window, &feature_config, cores)?;

    let clusterer = FeatureClusterer {
        tolerance: alignment_config.mz_tolerance.clone(),
        rt_tolerance: alignment_config.rt_tolerance,
    };

    let mut slots = prepare_slots(
        clusterer.cluster(collect_tagged(&mut datasets)),
        datasets.len(),
        alignment_config.rt_tolerance,
    );

    fill_all_missing(
        &mut slots,
        &datasets,
        alignment_config.eic_options,
        alignment_config.peak_options,
        cores,
    )?;

    let sample_names: Vec<String> = datasets.iter().map(|(name, _, _)| name.clone()).collect();
    let results = build_results(
        slots,
        alignment_config.min_samples,
        datasets.len(),
        &sample_names,
    );

    Ok(dedup(
        results,
        &alignment_config.mz_tolerance,
        alignment_config.rt_tolerance,
    ))
}

fn collect_tagged(datasets: &mut [SampleDataset]) -> Vec<TaggedFeature> {
    datasets
        .iter_mut()
        .enumerate()
        .flat_map(|(idx, (_, _, features))| {
            std::mem::take(features)
                .into_iter()
                .map(move |feature| TaggedFeature {
                    sample_index: idx,
                    feature,
                })
        })
        .collect()
}

fn prepare_slots(clusters: Vec<Cluster>, n_samples: usize, rt_tol: f64) -> Vec<ClusterSlot> {
    clusters
        .into_iter()
        .filter_map(|cluster| {
            let slots = assign_best_per_sample(cluster, n_samples);
            let bounds = compute_search_bounds(&slots, rt_tol)?;
            Some((slots, bounds))
        })
        .collect()
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
fn fill_all_missing(
    slots: &mut [ClusterSlot],
    datasets: &[SampleDataset],
    eic_options: EicOptions,
    peak_options: Option<FindPeaksOptions>,
    cores: usize,
) -> Result<(), AlignmentError> {
    let shared_slots: &[ClusterSlot] = slots;

    let per_sample: Vec<(usize, Vec<(usize, Feature)>)> = run_with_cores(cores, || {
        datasets
            .par_iter()
            .enumerate()
            .map(|(sample_idx, (_, source, _))| {
                let filled = match source {
                    SampleSourceKind::Mzml(path) => {
                        let mut mzml = open_mzml(path).map_err(|e| AlignmentError::Parse {
                            path: path.to_string_lossy().to_string(),
                            source: e.to_string(),
                        })?;
                        let mut reader = EicReader::Mzml(&mut mzml);
                        fill_sample(
                            shared_slots,
                            sample_idx,
                            &mut reader,
                            eic_options,
                            peak_options.clone(),
                        )?
                    }
                    SampleSourceKind::Ion(path) => {
                        let mut owned = open_ion(path).map_err(|e| AlignmentError::Parse {
                            path: path.to_string_lossy().to_string(),
                            source: e.to_string(),
                        })?;
                        let mut reader = EicReader::Ion(&mut owned);
                        fill_sample(
                            shared_slots,
                            sample_idx,
                            &mut reader,
                            eic_options,
                            peak_options.clone(),
                        )?
                    }
                };
                Ok((sample_idx, filled))
            })
            .collect::<Result<Vec<_>, AlignmentError>>()
    })?;

    for (sample_idx, filled) in per_sample {
        for (ci, feature) in filled {
            slots[ci].0[sample_idx] = Some(feature);
        }
    }

    Ok(())
}

fn fill_sample(
    slots: &[ClusterSlot],
    sample_idx: usize,
    reader: &mut EicReader,
    eic_options: EicOptions,
    peak_options: Option<FindPeaksOptions>,
) -> Result<Vec<(usize, Feature)>, AlignmentError> {
    let missing: Vec<usize> = slots
        .iter()
        .enumerate()
        .filter_map(|(ci, (s, _))| {
            if s[sample_idx].is_none() {
                Some(ci)
            } else {
                None
            }
        })
        .collect();

    if missing.is_empty() {
        return Ok(Vec::new());
    }

    let rt_min = missing
        .iter()
        .map(|&ci| slots[ci].1.rt_from)
        .fold(f64::INFINITY, f64::min);
    let rt_max = missing
        .iter()
        .map(|&ci| slots[ci].1.rt_to)
        .fold(f64::NEG_INFINITY, f64::max);

    let scan_times = get_scan_times(reader, rt_min, rt_max, MS1_LEVEL);

    if scan_times.is_empty() {
        return Ok(Vec::new());
    }

    let all_times: Vec<f64> = scan_times.iter().map(|s| s.rt).collect();

    let mut filled: Vec<(usize, Feature)> = Vec::new();

    let mut jobs = missing;
    jobs.sort_unstable_by(|&a, &b| {
        slots[a]
            .1
            .target_mz
            .partial_cmp(&slots[b].1.target_mz)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut job_idx = 0;
    while job_idx < jobs.len() {
        let tile_center_mz = slots[jobs[job_idx]].1.target_mz;
        let tile_tolerance = mz_tolerance_for(tile_center_mz, eic_options);

        let mut tile_jobs = vec![jobs[job_idx]];
        job_idx += 1;
        while job_idx < jobs.len()
            && (slots[jobs[job_idx]].1.target_mz - tile_center_mz).abs() <= tile_tolerance * 2.0
        {
            tile_jobs.push(jobs[job_idx]);
            job_idx += 1;
        }

        let tile_rt_min = tile_jobs
            .iter()
            .map(|&ci| slots[ci].1.rt_from)
            .fold(f64::INFINITY, f64::min);
        let tile_rt_max = tile_jobs
            .iter()
            .map(|&ci| slots[ci].1.rt_to)
            .fold(f64::NEG_INFINITY, f64::max);

        let start = lower_bound(&all_times, tile_rt_min);
        let end = upper_bound(&all_times, tile_rt_max);

        if start >= end {
            continue;
        }

        let tile_mz_lo = tile_center_mz - tile_tolerance;
        let tile_mz_hi = tile_jobs
            .iter()
            .map(|&ci| slots[ci].1.target_mz)
            .fold(f64::NEG_INFINITY, f64::max)
            + tile_tolerance;

        let mut tile_scans: Vec<(f64, Vec<f64>, Vec<f64>)> = Vec::new();
        for scan_time in scan_times[start..end].iter() {
            let mut mz = Vec::new();
            let mut intensity = Vec::new();
            read_mz_window(
                reader,
                scan_time.index,
                tile_mz_lo,
                tile_mz_hi,
                &mut mz,
                &mut intensity,
            )
            .map_err(|e| AlignmentError::FastPath {
                path: format!("scan {}", scan_time.index),
                source: e,
            })?;
            tile_scans.push((scan_time.rt, mz, intensity));
        }

        #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
        let tile_iter = tile_jobs.par_iter();
        #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
        let tile_iter = tile_jobs.iter();

        let tile_features: Vec<(usize, Feature)> = tile_iter
            .filter_map(|&ci| {
                let bounds = &slots[ci].1;
                let tolerance = mz_tolerance_for(bounds.target_mz, eic_options);
                let window_lo = bounds.target_mz - tolerance;
                let window_hi = bounds.target_mz + tolerance;

                let mut intensities = Vec::with_capacity(tile_scans.len());
                for (_, mz, intensity) in &tile_scans {
                    let mut sum = 0.0;
                    for (m, i) in mz.iter().zip(intensity.iter()) {
                        if *m >= window_lo && *m <= window_hi {
                            sum += i;
                        }
                    }
                    intensities.push(sum);
                }

                let time_slice: Vec<f64> = tile_scans.iter().map(|(rt, _, _)| *rt).collect();

                get_peak(
                    &DataXY {
                        x: time_slice,
                        y: intensities,
                    },
                    &Roi {
                        rt: bounds.center_rt,
                        half_width: bounds.rt_to - bounds.rt_from,
                    },
                    peak_options.clone(),
                )
                .filter(|p| p.intensity > 0.0)
                .map(|p| {
                    let measured_mz =
                        calculate_weighted_mz(&tile_scans, window_lo, window_hi, p.from, p.to)
                            .unwrap_or(bounds.target_mz);
                    let feature = Feature {
                        mz: measured_mz,
                        rt: p.rt,
                        intensity: p.intensity,
                        from: p.from,
                        to: p.to,
                        n_points: p.n_points,
                        integral: p.integral,
                        noise: p.noise,
                    };
                    (ci, feature)
                })
            })
            .collect();

        filled.extend(tile_features);
    }

    Ok(filled)
}

fn calculate_weighted_mz(
    scans: &[(f64, Vec<f64>, Vec<f64>)],
    mz_lo: f64,
    mz_hi: f64,
    rt_from: f64,
    rt_to: f64,
) -> Option<f64> {
    let mut total_weighted_mz = 0.0;
    let mut total_intensity = 0.0;

    for (rt, mz, intensity) in scans {
        if *rt < rt_from || *rt > rt_to {
            continue;
        }
        for (m, i) in mz.iter().zip(intensity.iter()) {
            if *m < mz_lo || *m > mz_hi {
                continue;
            }
            total_weighted_mz += m * i;
            total_intensity += i;
        }
    }

    if total_intensity > 0.0 {
        Some(total_weighted_mz / total_intensity)
    } else {
        None
    }
}

pub fn weighted_centroid_mz(
    scans: &[CentroidScan],
    times: &[f64],
    target_mz: f64,
    options: EicOptions,
    rt_from: f64,
    rt_to: f64,
) -> Option<f64> {
    if !target_mz.is_finite() || target_mz <= 0.0 {
        return None;
    }
    let tol = mz_tolerance_for(target_mz, options);
    if !tol.is_finite() || tol <= 0.0 {
        return None;
    }
    let (lo, hi) = (target_mz - tol, target_mz + tol);
    let (mut wsum, mut isum) = (0.0f64, 0.0f64);
    for (scan, &rt) in scans.iter().zip(times.iter()) {
        if rt < rt_from || rt > rt_to {
            continue;
        }
        let start = lower_bound(&scan.mz, lo);
        let end = upper_bound(&scan.mz, hi);
        for i in start..end {
            wsum += scan.mz[i] * scan.intensity[i];
            isum += scan.intensity[i];
        }
    }
    (isum > 0.0).then(|| wsum / isum)
}

fn build_results(
    slots: Vec<ClusterSlot>,
    min_samples: usize,
    total_samples: usize,
    sample_names: &[String],
) -> Vec<ConsensusFeature> {
    slots
        .into_iter()
        .filter_map(|(s, bounds)| {
            let debug_missing = collect_missing_samples(&s, &bounds, sample_names);
            let hits = collect_filled_slots(s);
            require_minimum_frequency(hits, min_samples).map(|hits| {
                let mut feature = aggregate_into_consensus(hits, &bounds, total_samples);
                feature.debug_missing = debug_missing;
                feature
            })
        })
        .collect()
}

pub(crate) fn dedup(
    mut results: Vec<ConsensusFeature>,
    tolerance: &MzTolerance,
    rt_tol: f64,
) -> Vec<ConsensusFeature> {
    results.sort_unstable_by(|a, b| {
        b.n_samples.cmp(&a.n_samples).then_with(|| {
            b.intensity
                .partial_cmp(&a.intensity)
                .unwrap_or(Ordering::Equal)
        })
    });
    let mut kept: Vec<ConsensusFeature> = Vec::new();
    'outer: for f in results {
        for k in &kept {
            if tolerance.are_close_to_ref(f.mz, k.mz) && (f.rt - k.rt).abs() <= rt_tol {
                continue 'outer;
            }
        }
        kept.push(f);
    }
    kept
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
fn load_sample_files(directory: &str) -> Result<Vec<(String, SampleSourceKind)>, AlignmentError> {
    let entries: Vec<_> = fs::read_dir(directory)?
        .filter_map(|entry| {
            let path = entry.ok()?.path();
            matches!(path.extension()?.to_str()?, "mzML" | "ion").then_some(path)
        })
        .collect();

    if entries.is_empty() {
        return Err(AlignmentError::NoFiles);
    }

    entries
        .into_iter()
        .map(|path| {
            let file_name = path
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown")
                .to_string();
            let ext = path.extension().and_then(|s| s.to_str()).unwrap_or("");
            let source = match ext {
                "mzML" => SampleSourceKind::Mzml(path),
                "ion" => SampleSourceKind::Ion(path),
                other => return Err(AlignmentError::UnsupportedFormat(other.to_string())),
            };
            Ok((file_name, source))
        })
        .collect()
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
fn detect_features_per_sample(
    samples: Vec<(String, SampleSourceKind)>,
    time_window: FromTo,
    config: &FindFeaturesOptions,
    cores: usize,
) -> Result<Vec<SampleDataset>, AlignmentError> {
    let mut results = Vec::new();

    for (idx, (name, source)) in samples.into_iter().enumerate() {
        let start = Instant::now();
        let features = match &source {
            SampleSourceKind::Mzml(path) => {
                let mut mzml = open_mzml(path).map_err(|e| AlignmentError::Parse {
                    path: path.to_string_lossy().to_string(),
                    source: e.to_string(),
                })?;
                let mut reader = EicReader::Mzml(&mut mzml);
                find_features(&mut reader, time_window, Some(config.clone()), cores).map_err(
                    |e| AlignmentError::FeatureDetection {
                        path: path.to_string_lossy().to_string(),
                        source: e,
                    },
                )?
            }
            SampleSourceKind::Ion(path) => {
                let mut owned = open_ion(path).map_err(|e| AlignmentError::Parse {
                    path: path.to_string_lossy().to_string(),
                    source: e.to_string(),
                })?;
                let mut reader = EicReader::Ion(&mut owned);
                find_features(&mut reader, time_window, Some(config.clone()), cores).map_err(
                    |e| AlignmentError::FeatureDetection {
                        path: path.to_string_lossy().to_string(),
                        source: e,
                    },
                )?
            }
        };
        eprintln!(
            "[Sample {}] {} processed in {:.3}s",
            idx,
            name,
            start.elapsed().as_secs_f64()
        );
        results.push((name, source, features));
    }

    Ok(results)
}

pub(crate) fn assign_best_per_sample(cluster: Cluster, n_samples: usize) -> Vec<Option<Feature>> {
    let mut slots: Vec<Option<Feature>> = vec![None; n_samples];
    for tagged in cluster {
        let idx = tagged.sample_index;
        if slots[idx]
            .as_ref()
            .is_none_or(|f| tagged.feature.intensity > f.intensity)
        {
            slots[idx] = Some(tagged.feature);
        }
    }
    slots
}

pub(crate) fn compute_search_bounds(
    slots: &[Option<Feature>],
    rt_tol: f64,
) -> Option<SearchBounds> {
    let present: Vec<&Feature> = slots.iter().flatten().collect();
    if present.is_empty() {
        return None;
    }
    let mut mzs: Vec<f64> = present.iter().map(|f| f.mz).collect();
    let mut rts: Vec<f64> = present.iter().map(|f| f.rt).collect();
    mzs.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    rts.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let med_mz = mzs[mzs.len() / 2];
    let med_rt = rts[rts.len() / 2];
    Some(SearchBounds {
        target_mz: med_mz,
        center_rt: med_rt,
        rt_from: med_rt - rt_tol,
        rt_to: med_rt + rt_tol,
    })
}

pub(crate) fn collect_filled_slots(slots: Vec<Option<Feature>>) -> Vec<Feature> {
    slots.into_iter().flatten().collect()
}

pub(crate) fn require_minimum_frequency(
    hits: Vec<Feature>,
    threshold: usize,
) -> Option<Vec<Feature>> {
    if hits.len() >= threshold {
        Some(hits)
    } else {
        None
    }
}

pub(crate) fn aggregate_into_consensus(
    hits: Vec<Feature>,
    bounds: &SearchBounds,
    total_samples: usize,
) -> ConsensusFeature {
    let n = hits.len();

    let mut mz_values: Vec<f64> = hits.iter().map(|f| f.mz).collect();
    let mut rt_values: Vec<f64> = hits.iter().map(|f| f.rt).collect();
    let mut intensity_values: Vec<f64> = hits.iter().map(|f| f.intensity).collect();
    let mut integral_values: Vec<f64> = hits.iter().map(|f| f.integral).collect();

    ConsensusFeature {
        mz: median(&mut mz_values),
        rt: median(&mut rt_values),
        from: bounds.rt_from,
        to: bounds.rt_to,
        intensity: median(&mut intensity_values),
        integral: median(&mut integral_values),
        frequency: if total_samples > 0 {
            n as f64 / total_samples as f64
        } else {
            0.0
        },
        n_samples: n,
        debug_missing: Vec::new(),
    }
}

fn collect_missing_samples(
    slots: &[Option<Feature>],
    bounds: &SearchBounds,
    sample_names: &[String],
) -> Vec<MissingSample> {
    slots
        .iter()
        .enumerate()
        .filter(|(_, slot)| slot.is_none())
        .map(|(index, _)| MissingSample {
            sample: sample_names.get(index).cloned().unwrap_or_default(),
            target_mz: bounds.target_mz,
            center_rt: bounds.center_rt,
            rt_from: bounds.rt_from,
            rt_to: bounds.rt_to,
        })
        .collect()
}

pub(crate) fn median(values: &mut [f64]) -> f64 {
    values.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let mid = values.len() / 2;
    if values.len().is_multiple_of(2) {
        (values[mid - 1] + values[mid]) / 2.0
    } else {
        values[mid]
    }
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
fn open_mzml(path: &Path) -> Result<MzML, String> {
    let file = File::open(path).map_err(|e| format!("open {}: {}", path.display(), e))?;
    let mmap =
        unsafe { Mmap::map(&file) }.map_err(|e| format!("mmap {}: {}", path.display(), e))?;
    parse_mzml(&mmap[..]).map_err(|e| e.to_string())
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
#[inline]
fn open_ion(path: &Path) -> Result<IonReader, String> {
    IonReader::open_file(
        path,
        ReadOptions {
            max_cached_bytes: ION_CACHE_BYTES,
            ..Default::default()
        },
    )
    .map_err(|e| e.to_string())
}
