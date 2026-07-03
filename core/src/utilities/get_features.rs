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
        CentroidScan, EicOptions, EicReader, MS1_LEVEL, SpectrumKind, get_scan_times,
        get_spectrum_kind, highest_intensity_in_window, lower_bound, mz_tolerance_for,
        read_mz_window, upper_bound,
    },
    find_features::{Feature, FeatureError, FindFeaturesOptions, MzTolerance, find_features},
    find_masses::find_masses,
    find_peaks::FindPeaksOptions,
    get_peak::get_peak,
    math::median,
    mz_estimator::{
        MzEstimator, MzEstimatorKind, SampleMz, has_multiple_masses, make_estimator,
        same_mass_gap, split_by_gap,
    },
    structs::{DataXY, FromTo, Peak, Roi},
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
    pub mz_estimator: MzEstimatorKind,
    pub fill_window_minutes: f64,
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
            mz_estimator: MzEstimatorKind::default(),
            fill_window_minutes: 2.0,
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
        let by_mz = |a: &TaggedFeature, b: &TaggedFeature| {
            a.feature
                .mz
                .partial_cmp(&b.feature.mz)
                .unwrap_or(Ordering::Equal)
        };

        #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
        tagged.par_sort_unstable_by(by_mz);
        #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
        tagged.sort_unstable_by(by_mz);

        let mz_groups = self.group_by_mz(tagged);

        #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
        {
            mz_groups
                .into_par_iter()
                .flat_map_iter(|group| self.subdivide_by_rt(group.into_items()))
                .filter_map(|growing| self.keep_cluster(growing))
                .collect()
        }
        #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
        {
            mz_groups
                .into_iter()
                .flat_map(|group| self.subdivide_by_rt(group.into_items()))
                .filter_map(|growing| self.keep_cluster(growing))
                .collect()
        }
    }

    fn keep_cluster(&self, growing: MzRtCluster) -> Option<Cluster> {
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

#[derive(Clone, Copy)]
pub(crate) struct MassPeak {
    pub(crate) mz: f64,
    pub(crate) intensity: f64,
    pub(crate) integral: f64,
    pub(crate) rt: f64,
}

struct ClusterSlot {
    features: Vec<Option<Feature>>,
    apex_values: Vec<Option<SampleMz>>,
    masses: Vec<Option<Vec<MassPeak>>>,
    bounds: SearchBounds,
}

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

    let tagged = collect_tagged(&mut datasets);
    let clusters = run_with_cores(cores, || clusterer.cluster(tagged));
    let mut slots = prepare_slots(clusters, datasets.len(), alignment_config.rt_tolerance);

    let estimator = make_estimator(&alignment_config.mz_estimator);

    resolve_samples(
        &mut slots,
        &datasets,
        alignment_config.eic_options,
        alignment_config.peak_options,
        alignment_config.fill_window_minutes,
        cores,
    )?;

    let sample_names: Vec<String> = datasets.iter().map(|(name, _, _)| name.clone()).collect();
    let (single, split) = build_results(
        slots,
        alignment_config.min_samples,
        datasets.len(),
        &sample_names,
        estimator.as_ref(),
        &alignment_config.mz_tolerance,
        alignment_config.eic_options.ppm_tolerance,
    );

    let mut results = run_with_cores(cores, || {
        dedup(
            single,
            &alignment_config.mz_tolerance,
            alignment_config.rt_tolerance,
        )
    });
    let same_mass = MzTolerance {
        mz_absolute: 0.0,
        ppm: 0.5 * alignment_config.eic_options.ppm_tolerance,
    };
    results.extend(run_with_cores(cores, || {
        dedup(split, &same_mass, alignment_config.rt_tolerance)
    }));
    Ok(results)
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
            let features = assign_best_per_sample(cluster, n_samples);
            let bounds = compute_search_bounds(&features, rt_tol)?;
            Some(ClusterSlot {
                apex_values: vec![None; features.len()],
                masses: vec![None; features.len()],
                features,
                bounds,
            })
        })
        .collect()
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
fn resolve_samples(
    slots: &mut [ClusterSlot],
    datasets: &[SampleDataset],
    eic_options: EicOptions,
    peak_options: Option<FindPeaksOptions>,
    fill_window: f64,
    cores: usize,
) -> Result<(), AlignmentError> {
    if slots.is_empty() {
        return Ok(());
    }

    let tiles = build_tiles(slots, eic_options);

    for (sample_idx, (_, source, _)) in datasets.iter().enumerate() {
        let loaded = match source {
            SampleSourceKind::Mzml(path) => {
                let mut mzml = open_mzml(path).map_err(|e| AlignmentError::Parse {
                    path: path.to_string_lossy().to_string(),
                    source: e.to_string(),
                })?;
                let mut reader = EicReader::Mzml(&mut mzml);
                load_sample_scans(&mut reader, slots, sample_idx, eic_options, fill_window)?
            }
            SampleSourceKind::Ion(path) => {
                let mut owned = open_ion(path).map_err(|e| AlignmentError::Parse {
                    path: path.to_string_lossy().to_string(),
                    source: e.to_string(),
                })?;
                let mut reader = EicReader::Ion(&mut owned);
                load_sample_scans(&mut reader, slots, sample_idx, eic_options, fill_window)?
            }
        };

        let Some((kind, all_times, scans)) = loaded else {
            continue;
        };

        let shared_slots: &[ClusterSlot] = slots;
        let items: Vec<(usize, Resolved)> = run_with_cores(cores, || {
            #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
            {
                tiles
                    .par_iter()
                    .flat_map_iter(|tile| {
                        resolve_tile(
                            tile,
                            shared_slots,
                            sample_idx,
                            &all_times,
                            &scans,
                            eic_options,
                            &peak_options,
                            fill_window,
                            kind,
                        )
                    })
                    .collect()
            }
            #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
            {
                tiles
                    .iter()
                    .flat_map(|tile| {
                        resolve_tile(
                            tile,
                            shared_slots,
                            sample_idx,
                            &all_times,
                            &scans,
                            eic_options,
                            &peak_options,
                            fill_window,
                            kind,
                        )
                    })
                    .collect()
            }
        });

        for (ci, resolved) in items {
            if let Some(feature) = resolved.feature {
                slots[ci].features[sample_idx] = Some(feature);
            }
            if let Some(apex) = resolved.apex {
                slots[ci].apex_values[sample_idx] = Some(apex);
            }
            if let Some(masses) = resolved.masses {
                slots[ci].masses[sample_idx] = Some(masses);
            }
        }
    }

    Ok(())
}

pub(crate) struct Resolved {
    pub(crate) feature: Option<Feature>,
    pub(crate) apex: Option<SampleMz>,
    pub(crate) masses: Option<Vec<MassPeak>>,
}

fn build_tiles(slots: &[ClusterSlot], eic_options: EicOptions) -> Vec<Vec<usize>> {
    let mut jobs: Vec<usize> = (0..slots.len()).collect();
    jobs.sort_unstable_by(|&a, &b| {
        slots[a]
            .bounds
            .target_mz
            .partial_cmp(&slots[b].bounds.target_mz)
            .unwrap_or(Ordering::Equal)
    });

    let mut tiles: Vec<Vec<usize>> = Vec::new();
    let mut job_idx = 0;
    while job_idx < jobs.len() {
        let tile_center_mz = slots[jobs[job_idx]].bounds.target_mz;
        let tile_tolerance = mz_tolerance_for(tile_center_mz, eic_options);

        let mut tile = vec![jobs[job_idx]];
        job_idx += 1;
        while job_idx < jobs.len()
            && (slots[jobs[job_idx]].bounds.target_mz - tile_center_mz).abs() <= tile_tolerance * 2.0
        {
            tile.push(jobs[job_idx]);
            job_idx += 1;
        }
        tiles.push(tile);
    }
    tiles
}

type LoadedScans = Option<(SpectrumKind, Vec<f64>, Vec<(f64, Vec<f64>, Vec<f64>)>)>;

fn load_sample_scans(
    reader: &mut EicReader,
    slots: &[ClusterSlot],
    sample_idx: usize,
    eic_options: EicOptions,
    fill_window: f64,
) -> Result<LoadedScans, AlignmentError> {
    let kind = get_spectrum_kind(reader);
    let rt_min = slots
        .iter()
        .map(|slot| sample_rt_from(slot, sample_idx, fill_window))
        .fold(f64::INFINITY, f64::min);
    let rt_max = slots
        .iter()
        .map(|slot| sample_rt_to(slot, sample_idx, fill_window))
        .fold(f64::NEG_INFINITY, f64::max);

    let scan_times = get_scan_times(reader, rt_min, rt_max, MS1_LEVEL);
    if scan_times.is_empty() {
        return Ok(None);
    }

    let mut mz_lo_all = f64::INFINITY;
    let mut mz_hi_all = f64::NEG_INFINITY;
    for slot in slots {
        let tolerance = mz_tolerance_for(slot.bounds.target_mz, eic_options);
        mz_lo_all = mz_lo_all.min(slot.bounds.target_mz - tolerance);
        mz_hi_all = mz_hi_all.max(slot.bounds.target_mz + tolerance);
    }

    let mut scans: Vec<(f64, Vec<f64>, Vec<f64>)> = Vec::with_capacity(scan_times.len());
    let mut mz = Vec::new();
    let mut intensity = Vec::new();
    for scan_time in &scan_times {
        read_mz_window(
            reader,
            scan_time.index,
            mz_lo_all,
            mz_hi_all,
            &mut mz,
            &mut intensity,
        )
        .map_err(|e| AlignmentError::FastPath {
            path: format!("scan {}", scan_time.index),
            source: e,
        })?;
        scans.push((scan_time.rt, mz.clone(), intensity.clone()));
    }

    let all_times: Vec<f64> = scans.iter().map(|s| s.0).collect();
    Ok(Some((kind, all_times, scans)))
}

#[allow(clippy::too_many_arguments)]
fn resolve_tile(
    tile: &[usize],
    slots: &[ClusterSlot],
    sample_idx: usize,
    all_times: &[f64],
    scans: &[(f64, Vec<f64>, Vec<f64>)],
    eic_options: EicOptions,
    peak_options: &Option<FindPeaksOptions>,
    fill_window: f64,
    kind: SpectrumKind,
) -> Vec<(usize, Resolved)> {
    let tile_center_mz = slots[tile[0]].bounds.target_mz;
    let tile_tolerance = mz_tolerance_for(tile_center_mz, eic_options);

    let tile_rt_min = tile
        .iter()
        .map(|&ci| sample_rt_from(&slots[ci], sample_idx, fill_window))
        .fold(f64::INFINITY, f64::min);
    let tile_rt_max = tile
        .iter()
        .map(|&ci| sample_rt_to(&slots[ci], sample_idx, fill_window))
        .fold(f64::NEG_INFINITY, f64::max);

    let start = lower_bound(all_times, tile_rt_min);
    let end = upper_bound(all_times, tile_rt_max);
    if start >= end {
        return Vec::new();
    }

    let tile_mz_lo = tile_center_mz - tile_tolerance;
    let tile_mz_hi = tile
        .iter()
        .map(|&ci| slots[ci].bounds.target_mz)
        .fold(f64::NEG_INFINITY, f64::max)
        + tile_tolerance;

    let tile_scans: Vec<(f64, Vec<f64>, Vec<f64>)> = scans[start..end]
        .iter()
        .map(|(rt, mz, intensity)| {
            let lo = lower_bound(mz, tile_mz_lo);
            let hi = upper_bound(mz, tile_mz_hi).min(mz.len());
            (*rt, mz[lo..hi].to_vec(), intensity[lo..hi].to_vec())
        })
        .collect();

    tile.iter()
        .map(|&ci| {
            let bounds = &slots[ci].bounds;
            let tolerance = mz_tolerance_for(bounds.target_mz, eic_options);
            let mz_lo = bounds.target_mz - tolerance;
            let mz_hi = bounds.target_mz + tolerance;
            let existing = slots[ci].features[sample_idx].as_ref();
            let resolved = resolve_cluster(
                &tile_scans,
                existing,
                bounds,
                mz_lo,
                mz_hi,
                eic_options.ppm_tolerance,
                peak_options.clone(),
                kind,
            );
            (ci, resolved)
        })
        .collect()
}

fn sample_rt_from(slot: &ClusterSlot, sample_idx: usize, fill_window: f64) -> f64 {
    match slot.features[sample_idx].as_ref() {
        Some(feature) => feature.from,
        None => (slot.bounds.center_rt - fill_window).min(slot.bounds.rt_from),
    }
}

fn sample_rt_to(slot: &ClusterSlot, sample_idx: usize, fill_window: f64) -> f64 {
    match slot.features[sample_idx].as_ref() {
        Some(feature) => feature.to,
        None => (slot.bounds.center_rt + fill_window).max(slot.bounds.rt_to),
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn resolve_cluster(
    tile_scans: &[(f64, Vec<f64>, Vec<f64>)],
    existing: Option<&Feature>,
    bounds: &SearchBounds,
    mz_lo: f64,
    mz_hi: f64,
    ppm: f64,
    peak_options: Option<FindPeaksOptions>,
    kind: SpectrumKind,
) -> Resolved {
    let (from, to, filled) = match existing {
        Some(feature) => (feature.from, feature.to, None),
        None => match detect_peak(tile_scans, bounds, mz_lo, mz_hi, peak_options.clone()) {
            Some(feature) => (feature.from, feature.to, Some(feature)),
            None => {
                return Resolved {
                    feature: None,
                    apex: None,
                    masses: None,
                };
            }
        },
    };

    let window = crop_scans(tile_scans, from, to);
    let mass_points = find_apex_masses(window, mz_lo, mz_hi, kind);
    let apex = mass_points
        .iter()
        .copied()
        .max_by(|a, b| a.intensity.partial_cmp(&b.intensity).unwrap_or(Ordering::Equal));
    let peak_feature = existing.or(filled.as_ref());

    let masses: Vec<MassPeak> = if mass_points.len() <= 1 {
        match (mass_points.first(), peak_feature) {
            (Some(centroid), Some(feature)) => vec![MassPeak {
                mz: centroid.mz,
                intensity: feature.intensity,
                integral: feature.integral,
                rt: feature.rt,
            }],
            _ => Vec::new(),
        }
    } else {
        mass_points
            .iter()
            .filter_map(|centroid| {
                let gap = same_mass_gap(centroid.mz, ppm);
                measure_mass(window, centroid.mz, gap, bounds, peak_options.clone()).map(|peak| {
                    MassPeak {
                        mz: centroid.mz,
                        intensity: peak.intensity,
                        integral: peak.integral,
                        rt: peak.rt,
                    }
                })
            })
            .collect()
    };

    let feature = filled.map(|mut feature| {
        if let Some(sample_mz) = apex {
            feature.mz = sample_mz.mz;
        }
        feature
    });

    Resolved {
        feature,
        apex,
        masses: (!masses.is_empty()).then_some(masses),
    }
}

fn measure_mass(
    window: &[(f64, Vec<f64>, Vec<f64>)],
    center_mz: f64,
    gap: f64,
    bounds: &SearchBounds,
    peak_options: Option<FindPeaksOptions>,
) -> Option<Peak> {
    let mz_lo = center_mz - gap;
    let mz_hi = center_mz + gap;
    let mut times = Vec::with_capacity(window.len());
    let mut intensities = Vec::with_capacity(window.len());
    for (rt, mz, intensity) in window {
        times.push(*rt);
        intensities.push(highest_intensity_in_window(mz, intensity, mz_lo, mz_hi));
    }
    anchor_edges(&mut intensities);
    get_peak(
        &DataXY {
            x: times,
            y: intensities,
        },
        &Roi {
            rt: bounds.center_rt,
            half_width: bounds.rt_to - bounds.rt_from,
        },
        peak_options,
    )
    .filter(|p| p.intensity > 0.0)
}

fn find_apex_masses(
    scans: &[(f64, Vec<f64>, Vec<f64>)],
    mz_lo: f64,
    mz_hi: f64,
    kind: SpectrumKind,
) -> Vec<SampleMz> {
    let mut masses: Vec<SampleMz> = Vec::new();
    let mut best_intensity = 0.0;
    for (_, mz, intensity) in scans {
        let top = highest_intensity_in_window(mz, intensity, mz_lo, mz_hi);
        if top > best_intensity {
            best_intensity = top;
            masses = find_masses(mz, intensity, mz_lo, mz_hi, kind);
        }
    }
    masses
}

fn detect_peak(
    tile_scans: &[(f64, Vec<f64>, Vec<f64>)],
    bounds: &SearchBounds,
    mz_lo: f64,
    mz_hi: f64,
    peak_options: Option<FindPeaksOptions>,
) -> Option<Feature> {
    let mut times = Vec::with_capacity(tile_scans.len());
    let mut intensities = Vec::with_capacity(tile_scans.len());
    for (rt, mz, intensity) in tile_scans {
        times.push(*rt);
        intensities.push(highest_intensity_in_window(mz, intensity, mz_lo, mz_hi));
    }

    anchor_edges(&mut intensities);

    get_peak(
        &DataXY {
            x: times,
            y: intensities,
        },
        &Roi {
            rt: bounds.center_rt,
            half_width: bounds.rt_to - bounds.rt_from,
        },
        peak_options,
    )
    .filter(|p| p.intensity > 0.0)
    .map(|p| Feature {
        mz: bounds.target_mz,
        rt: p.rt,
        intensity: p.intensity,
        from: p.from,
        to: p.to,
        n_points: p.n_points,
        integral: p.integral,
        noise: p.noise,
    })
}

const EDGE_BLEED_FRACTION: f64 = 0.05;

fn anchor_edges(values: &mut [f64]) {
    let length = values.len();
    if length < 3 {
        return;
    }
    let baseline = lowest_value(values);
    let limit = length / 3;
    flatten_leading_bleed(values, baseline, limit);
    flatten_trailing_bleed(values, baseline, limit);
}

fn lowest_value(values: &[f64]) -> f64 {
    let mut lowest = f64::INFINITY;
    for &value in values {
        if value.is_finite() && value < lowest {
            lowest = value;
        }
    }
    if lowest.is_finite() { lowest } else { 0.0 }
}

fn flatten_leading_bleed(values: &mut [f64], baseline: f64, limit: usize) {
    let threshold = baseline + EDGE_BLEED_FRACTION * (values[0] - baseline);
    let mut index = 0;
    while index < limit && values[index] > threshold {
        index += 1;
    }
    if index < limit {
        for value in values.iter_mut().take(index) {
            *value = baseline;
        }
    }
}

fn flatten_trailing_bleed(values: &mut [f64], baseline: f64, limit: usize) {
    let length = values.len();
    let threshold = baseline + EDGE_BLEED_FRACTION * (values[length - 1] - baseline);
    let mut index = length;
    while index > length - limit && values[index - 1] > threshold {
        index -= 1;
    }
    if index > length - limit {
        for value in values.iter_mut().skip(index) {
            *value = baseline;
        }
    }
}

fn crop_scans(
    scans: &[(f64, Vec<f64>, Vec<f64>)],
    from: f64,
    to: f64,
) -> &[(f64, Vec<f64>, Vec<f64>)] {
    let start = scans.partition_point(|(rt, _, _)| *rt < from);
    let end = scans.partition_point(|(rt, _, _)| *rt <= to);
    &scans[start..end.max(start)]
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

struct BuildContext<'a> {
    estimator: &'a dyn MzEstimator,
    tolerance: &'a MzTolerance,
    total_samples: usize,
}

#[allow(clippy::too_many_arguments)]
fn build_results(
    slots: Vec<ClusterSlot>,
    min_samples: usize,
    total_samples: usize,
    sample_names: &[String],
    estimator: &(dyn MzEstimator + Send + Sync),
    tolerance: &MzTolerance,
    ppm: f64,
) -> (Vec<ConsensusFeature>, Vec<ConsensusFeature>) {
    let context = BuildContext {
        estimator,
        tolerance,
        total_samples,
    };
    let mut single: Vec<ConsensusFeature> = Vec::new();
    let mut split: Vec<ConsensusFeature> = Vec::new();
    for slot in slots {
        let ClusterSlot {
            features,
            apex_values,
            masses,
            bounds,
        } = slot;
        let debug_missing = collect_missing_samples(&features, &bounds, sample_names);

        let mut pool: Vec<(usize, MassPeak)> = Vec::new();
        for (sample_idx, observations) in masses.iter().enumerate() {
            if let Some(observations) = observations {
                for observation in observations {
                    pool.push((sample_idx, *observation));
                }
            }
        }
        let pool_mz: Vec<f64> = pool.iter().map(|(_, o)| o.mz).collect();
        let cutoff = same_mass_gap(bounds.target_mz, ppm);

        if pool.is_empty() || !has_multiple_masses(&pool_mz, cutoff) {
            let dominant_apex: Vec<SampleMz> = apex_values.iter().flatten().copied().collect();
            let hits = collect_filled_slots(features);
            if let Some(hits) = require_minimum_frequency(hits, min_samples) {
                let mut feature = aggregate_into_consensus(hits, &bounds, total_samples);
                if let Some(mz) = context.estimator.combine(&dominant_apex, context.tolerance) {
                    feature.mz = mz;
                }
                feature.debug_missing = debug_missing;
                single.push(feature);
            }
            continue;
        }

        let mut group_features: Vec<ConsensusFeature> = Vec::new();
        for group in split_by_gap(&pool_mz, cutoff) {
            let low = group.iter().copied().fold(f64::INFINITY, f64::min);
            let high = group.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let entries: Vec<(usize, MassPeak)> = pool
                .iter()
                .copied()
                .filter(|(_, o)| o.mz >= low && o.mz <= high)
                .collect();
            if let Some(feature) = build_mass_feature(&entries, &bounds, &context)
                && feature.n_samples >= min_samples
            {
                group_features.push(feature);
            }
        }
        if let Some(first) = group_features.first_mut() {
            first.debug_missing = debug_missing;
        }
        split.extend(group_features);
    }
    (single, split)
}

fn build_mass_feature(
    entries: &[(usize, MassPeak)],
    bounds: &SearchBounds,
    context: &BuildContext,
) -> Option<ConsensusFeature> {
    if entries.is_empty() {
        return None;
    }

    let mut samples: Vec<usize> = entries.iter().map(|(sample_idx, _)| *sample_idx).collect();
    samples.sort_unstable();
    samples.dedup();

    let mut mzs: Vec<f64> = entries.iter().map(|(_, o)| o.mz).collect();
    let mut intensities: Vec<f64> = entries.iter().map(|(_, o)| o.intensity).collect();
    let mut integrals: Vec<f64> = entries.iter().map(|(_, o)| o.integral).collect();
    let mut rts: Vec<f64> = entries.iter().map(|(_, o)| o.rt).collect();

    Some(ConsensusFeature {
        mz: median(&mut mzs),
        rt: median(&mut rts),
        from: bounds.rt_from,
        to: bounds.rt_to,
        intensity: median(&mut intensities),
        integral: median(&mut integrals),
        frequency: if context.total_samples > 0 {
            samples.len() as f64 / context.total_samples as f64
        } else {
            0.0
        },
        n_samples: samples.len(),
        debug_missing: Vec::new(),
    })
}

pub(crate) fn dedup(
    results: Vec<ConsensusFeature>,
    tolerance: &MzTolerance,
    rt_tol: f64,
) -> Vec<ConsensusFeature> {
    if results.is_empty() {
        return results;
    }

    let mut by_mz = results;
    let by_mz_cmp =
        |a: &ConsensusFeature, b: &ConsensusFeature| a.mz.partial_cmp(&b.mz).unwrap_or(Ordering::Equal);
    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    by_mz.par_sort_unstable_by(by_mz_cmp);
    #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
    by_mz.sort_unstable_by(by_mz_cmp);

    let bands = split_dedup_bands(by_mz, tolerance);

    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    let mut kept: Vec<ConsensusFeature> = bands
        .into_par_iter()
        .flat_map_iter(|band| dedup_band(band, tolerance, rt_tol))
        .collect();
    #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
    let mut kept: Vec<ConsensusFeature> = bands
        .into_iter()
        .flat_map(|band| dedup_band(band, tolerance, rt_tol))
        .collect();

    #[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
    kept.par_sort_unstable_by(dedup_priority);
    #[cfg(all(target_arch = "wasm32", not(target_os = "wasi")))]
    kept.sort_unstable_by(dedup_priority);
    kept
}

fn dedup_priority(a: &ConsensusFeature, b: &ConsensusFeature) -> Ordering {
    b.n_samples.cmp(&a.n_samples).then_with(|| {
        b.intensity
            .partial_cmp(&a.intensity)
            .unwrap_or(Ordering::Equal)
    })
}

fn split_dedup_bands(
    sorted_by_mz: Vec<ConsensusFeature>,
    tolerance: &MzTolerance,
) -> Vec<Vec<ConsensusFeature>> {
    let mut bands: Vec<Vec<ConsensusFeature>> = Vec::new();
    let mut current: Vec<ConsensusFeature> = Vec::new();
    for feature in sorted_by_mz {
        if let Some(previous) = current.last() {
            let independent = !tolerance.are_close_to_ref(feature.mz, previous.mz)
                && !tolerance.are_close_to_ref(previous.mz, feature.mz);
            if independent {
                bands.push(std::mem::take(&mut current));
            }
        }
        current.push(feature);
    }
    if !current.is_empty() {
        bands.push(current);
    }
    bands
}

fn dedup_band(
    mut band: Vec<ConsensusFeature>,
    tolerance: &MzTolerance,
    rt_tol: f64,
) -> Vec<ConsensusFeature> {
    band.sort_unstable_by(dedup_priority);
    let mut kept: Vec<ConsensusFeature> = Vec::new();
    'outer: for feature in band {
        for existing in &kept {
            if tolerance.are_close_to_ref(feature.mz, existing.mz)
                && (feature.rt - existing.rt).abs() <= rt_tol
            {
                continue 'outer;
            }
        }
        kept.push(feature);
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

#[cfg(test)]
mod anchor_edges_tests {
    use super::anchor_edges;

    #[test]
    fn flattens_leading_bleed_and_keeps_center() {
        let mut values = vec![0.0; 30];
        let descent = [150_000.0, 90_000.0, 40_000.0, 10_000.0, 1_000.0];
        for (index, value) in descent.iter().enumerate() {
            values[index] = *value;
        }
        values[15] = 350_000.0;
        anchor_edges(&mut values);
        assert_eq!(values[0], 0.0);
        assert_eq!(values[1], 0.0);
        assert_eq!(values[15], 350_000.0);
    }

    #[test]
    fn flattens_trailing_bleed_and_keeps_center() {
        let mut values = vec![0.0; 30];
        let ascent = [1_000.0, 10_000.0, 40_000.0, 90_000.0, 150_000.0];
        for (offset, value) in ascent.iter().enumerate() {
            values[25 + offset] = *value;
        }
        values[15] = 350_000.0;
        anchor_edges(&mut values);
        assert_eq!(values[29], 0.0);
        assert_eq!(values[28], 0.0);
        assert_eq!(values[15], 350_000.0);
    }

    #[test]
    fn leaves_clean_edges_untouched() {
        let mut values = vec![0.0; 30];
        values[14] = 200_000.0;
        values[15] = 350_000.0;
        values[16] = 200_000.0;
        let before = values.clone();
        anchor_edges(&mut values);
        assert_eq!(values, before);
    }
}
