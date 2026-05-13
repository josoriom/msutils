use std::{
    cmp::Ordering,
    fmt::{Display, Formatter},
    io::Error,
    sync::Arc,
};

use rayon::prelude::*;

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use crate::utilities::gpu::GpuContext;

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use ionic::{
    ScanSource,
    ion::{DecoderConfig, Ion, OwnedIon},
};

use crate::utilities::{
    calculate_eic::{
        CentroidScan, EicOptions, ScanQuery, TimeUnit, get_eic_for_mz, get_scans, lower_bound,
        mz_tolerance_for, upper_bound,
    },
    find_features::{Feature, FeatureError, FindFeaturesOptions, MzTolerance, find_features},
    find_peaks::FindPeaksOptions,
    get_peak::get_peak,
    structs::{DataXY, FromTo, Roi},
};

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
const ION_CACHE_BYTES: usize = 128 * 1024 * 1024;

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use ionic::{MzML, parse_mzml};
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
    Parse { path: String, source: String },
    UnsupportedFormat(String),
    NoFiles,
    FeatureDetection(FeatureError),
}

impl Display for AlignmentError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Io(e) => write!(f, "IO error: {}", e),
            Self::Parse { path, source } => write!(f, "parse error for {}: {}", path, source),
            Self::UnsupportedFormat(s) => write!(f, "unsupported file format: {}", s),
            Self::NoFiles => write!(f, "no valid files found in directory"),
            Self::FeatureDetection(e) => write!(f, "feature detection error: {}", e),
        }
    }
}

impl std::error::Error for AlignmentError {}

impl From<Error> for AlignmentError {
    fn from(e: Error) -> Self {
        Self::Io(e)
    }
}

impl From<FeatureError> for AlignmentError {
    fn from(e: FeatureError) -> Self {
        Self::FeatureDetection(e)
    }
}

#[derive(Clone, Debug)]
pub struct ConsensusFeature {
    pub mz: f64,
    pub rmz: f64,
    pub rint: f64,
    pub rt: f64,
    pub intensity: f64,
    pub rintensity: f64,
    pub from: f64,
    pub to: f64,
    pub np: usize,
    pub integral: f64,
    pub frequency: usize,
}

#[derive(Clone, Debug)]
pub struct ConsensusAlignmentConfig {
    pub tolerance: MzTolerance,
    pub rt_tolerance: f64,
    pub frequency: usize,
    pub eic_options: EicOptions,
    pub peak_options: Option<FindPeaksOptions>,
}

impl Default for ConsensusAlignmentConfig {
    fn default() -> Self {
        Self {
            tolerance: MzTolerance {
                mz_abs: 0.005,
                ppm: 20.0,
            },
            rt_tolerance: 0.05,
            frequency: 1,
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
    pub(crate) seed_rt: f64,
    pub(crate) rt_from: f64,
    pub(crate) rt_to: f64,
}

#[derive(Debug)]
pub(crate) struct GrowingCluster {
    items: Vec<TaggedFeature>,
    sorted_mzs: Vec<f64>,
    sorted_rts: Vec<f64>,
    pub(crate) cached_median_mz: f64,
    pub(crate) cached_median_rt: f64,
}

impl GrowingCluster {
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

    fn group_by_mz(&self, items: Vec<TaggedFeature>) -> Vec<GrowingCluster> {
        let mut groups: Vec<GrowingCluster> = Vec::new();
        for item in items {
            let belongs = groups.last().is_some_and(|g| {
                self.tolerance
                    .are_close_to_ref(item.feature.mz, g.cached_median_mz)
            });
            if belongs {
                groups.last_mut().unwrap().push(item);
            } else {
                groups.push(GrowingCluster::new(item));
            }
        }
        groups
    }

    fn subdivide_by_rt(&self, group: Vec<TaggedFeature>) -> Vec<GrowingCluster> {
        let mut rt_sorted = group;
        rt_sorted.sort_unstable_by(|a, b| {
            a.feature
                .rt
                .partial_cmp(&b.feature.rt)
                .unwrap_or(Ordering::Equal)
        });
        let mut clusters: Vec<GrowingCluster> = Vec::new();
        for item in rt_sorted {
            let belongs = clusters
                .last()
                .is_some_and(|g| (item.feature.rt - g.cached_median_rt).abs() <= self.rt_tolerance);
            if belongs {
                clusters.last_mut().unwrap().push(item);
            } else {
                clusters.push(GrowingCluster::new(item));
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
    alignment_config: ConsensusAlignmentConfig,
    cores: usize,
) -> Result<Vec<ConsensusFeature>, AlignmentError> {
    let datasets = load_sample_files(directory_path)?;

    let mut feature_config = feature_config;
    if feature_config.use_gpu && feature_config.gpu_context.is_none() {
        feature_config.gpu_context = GpuContext::try_init().map(Arc::new);
    }

    let mut datasets = detect_features_per_sample(datasets, time_window, &feature_config, cores);

    let clusterer = FeatureClusterer {
        tolerance: alignment_config.tolerance.clone(),
        rt_tolerance: alignment_config.rt_tolerance,
    };

    let mut slots = prepare_slots(
        clusterer.cluster(collect_tagged(&mut datasets)),
        datasets.len(),
        alignment_config.rt_tolerance,
    );

    fill_all_missing(
        &mut slots,
        &mut datasets,
        alignment_config.eic_options,
        alignment_config.peak_options,
    );

    let results = build_results(slots, alignment_config.frequency);

    Ok(dedup(
        results,
        &alignment_config.tolerance,
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
    datasets: &mut [SampleDataset],
    eic_options: EicOptions,
    peak_options: Option<FindPeaksOptions>,
) {
    for (sample_idx, (name, source, _)) in datasets.iter_mut().enumerate() {
        match source {
            SampleSourceKind::Mzml(path) => match open_mzml(path) {
                Ok(mut mzml) => {
                    fill_sample(
                        slots,
                        sample_idx,
                        &mut mzml,
                        eic_options,
                        peak_options.clone(),
                    );
                }
                Err(e) => {
                    eprintln!("[Sample {}] {} failed to reopen: {}", sample_idx, name, e);
                }
            },
            SampleSourceKind::Ion(path) => match open_ion(path) {
                Ok(mut owned) => {
                    fill_sample(
                        slots,
                        sample_idx,
                        &mut *owned,
                        eic_options,
                        peak_options.clone(),
                    );
                }
                Err(e) => {
                    eprintln!("[Sample {}] {} failed to reopen: {}", sample_idx, name, e);
                }
            },
        }
    }
}

fn fill_sample(
    slots: &mut [ClusterSlot],
    sample_idx: usize,
    source: &mut impl ScanSource,
    eic_options: EicOptions,
    peak_options: Option<FindPeaksOptions>,
) {
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
        return;
    }

    let rt_min = missing
        .iter()
        .map(|&ci| slots[ci].1.rt_from)
        .fold(f64::INFINITY, f64::min);
    let rt_max = missing
        .iter()
        .map(|&ci| slots[ci].1.rt_to)
        .fold(f64::NEG_INFINITY, f64::max);

    let (all_times, all_scans) = get_scans(
        source,
        ScanQuery::RtRange(FromTo {
            from: rt_min,
            to: rt_max,
        }),
        TimeUnit::Minutes,
        1,
    );

    if all_scans.is_empty() {
        return;
    }

    let filled: Vec<(usize, Option<Feature>)> = missing
        .par_iter()
        .map(|&ci| {
            let bounds = &slots[ci].1;
            let start = lower_bound(&all_times, bounds.rt_from);
            let end = upper_bound(&all_times, bounds.rt_to);
            if start >= end {
                return (ci, None);
            }
            let time_slice = all_times[start..end].to_vec();
            let intensities = get_eic_for_mz(
                &all_scans[start..end],
                time_slice.len(),
                bounds.target_mz,
                eic_options,
            );
            let feature = get_peak(
                &DataXY {
                    x: time_slice,
                    y: intensities,
                },
                &Roi {
                    rt: bounds.seed_rt,
                    window: bounds.rt_to - bounds.rt_from,
                },
                peak_options.clone(),
            )
            .filter(|p| p.intensity > 0.0)
            .map(|p| {
                let measured_mz = weighted_centroid_mz(
                    &all_scans[start..end],
                    &all_times[start..end],
                    bounds.target_mz,
                    eic_options,
                    p.from,
                    p.to,
                )
                .unwrap_or(bounds.target_mz);
                Feature {
                    mz: measured_mz,
                    rt: p.rt,
                    intensity: p.intensity,
                    from: p.from,
                    to: p.to,
                    np: p.np,
                    integral: p.integral,
                    noise: p.noise,
                }
            });
            (ci, feature)
        })
        .collect();

    for (ci, f) in filled {
        if let Some(f) = f {
            slots[ci].0[sample_idx] = Some(f);
        }
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

fn build_results(slots: Vec<ClusterSlot>, frequency: usize) -> Vec<ConsensusFeature> {
    slots
        .into_iter()
        .filter_map(|(s, bounds)| {
            let hits = collect_filled_slots(s);
            require_minimum_frequency(hits, frequency)
                .map(|hits| aggregate_into_consensus(hits, &bounds))
        })
        .collect()
}

pub(crate) fn dedup(
    mut results: Vec<ConsensusFeature>,
    tolerance: &MzTolerance,
    rt_tol: f64,
) -> Vec<ConsensusFeature> {
    results.sort_unstable_by(|a, b| {
        b.frequency.cmp(&a.frequency).then_with(|| {
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
) -> Vec<SampleDataset> {
    samples
        .into_iter()
        .enumerate()
        .map(|(idx, (name, source))| {
            let start = Instant::now();
            let features = match &source {
                SampleSourceKind::Mzml(path) => match open_mzml(path) {
                    Ok(mut mzml) => {
                        find_features(&mut mzml, time_window, Some(config.clone()), cores)
                            .unwrap_or_default()
                    }
                    Err(e) => {
                        eprintln!("[Sample {}] {} failed to open: {}", idx, name, e);
                        Vec::new()
                    }
                },
                SampleSourceKind::Ion(path) => match open_ion(path) {
                    Ok(mut owned) => {
                        find_features(&mut *owned, time_window, Some(config.clone()), cores)
                            .unwrap_or_default()
                    }
                    Err(e) => {
                        eprintln!("[Sample {}] {} failed to open: {}", idx, name, e);
                        Vec::new()
                    }
                },
            };
            eprintln!(
                "[Sample {}] {} processed in {:.3}s",
                idx,
                name,
                start.elapsed().as_secs_f64()
            );
            (name, source, features)
        })
        .collect()
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
        seed_rt: med_rt,
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
) -> ConsensusFeature {
    let n = hits.len();

    let mut mz_values: Vec<f64> = hits.iter().map(|f| f.mz).collect();
    let mut rt_values: Vec<f64> = hits.iter().map(|f| f.rt).collect();
    let mut intensity_values: Vec<f64> = hits.iter().map(|f| f.intensity).collect();
    let mut np_values: Vec<f64> = hits.iter().map(|f| f.np as f64).collect();
    let mut integral_values: Vec<f64> = hits.iter().map(|f| f.integral).collect();

    ConsensusFeature {
        mz: median(&mut mz_values),
        rmz: rsd(&mz_values),
        rt: median(&mut rt_values),
        intensity: median(&mut intensity_values),
        rintensity: rsd(&intensity_values),
        from: bounds.rt_from,
        to: bounds.rt_to,
        np: median(&mut np_values) as usize,
        integral: median(&mut integral_values),
        rint: rsd(&integral_values),
        frequency: n,
    }
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

pub(crate) fn rsd(values: &[f64]) -> f64 {
    if values.len() < 2 {
        return 0.0;
    }
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    if mean == 0.0 {
        return 0.0;
    }
    let variance = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0);
    variance.sqrt() / mean
}

// TODO: Mmap-backed mzML loading is still being tested.
#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
fn open_mzml(path: &Path) -> Result<MzML, String> {
    let file = File::open(path).map_err(|e| format!("open {}: {}", path.display(), e))?;
    let mmap =
        unsafe { Mmap::map(&file) }.map_err(|e| format!("mmap {}: {}", path.display(), e))?;
    parse_mzml(&mmap[..]).map_err(|e| e.to_string())
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
#[inline]
fn open_ion(path: &Path) -> Result<OwnedIon, String> {
    Ion::open_file(
        path,
        DecoderConfig {
            max_cached_bytes: ION_CACHE_BYTES,
            ..Default::default()
        },
    )
    .map_err(|e| e.to_string())
}
