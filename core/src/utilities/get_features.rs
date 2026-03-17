use std::{
    cmp::Ordering,
    collections::HashSet,
    fmt::{Display, Formatter},
    io::Error,
};

use crate::utilities::find_peaks::FindPeaksOptions;
use crate::utilities::structs::FromTo;
use crate::utilities::{
    calculate_eic::{EicOptions, TimeUnit, collect_scans, compute_eic_for_mz},
    find_features::{
        Feature, FeatureError, FindFeaturesConfig, FindFeaturesOptions, MzTolerance, dedup_points,
        find_features,
    },
    get_peak::get_peak,
    structs::{DataXY, Roi},
};
use octo::MzML;

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use octo::{decoder::decode, parse_mzml};
#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
use std::{fs, time::Instant};

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
type SampleDataset = (String, MzML, Vec<Feature>);

pub(crate) struct SearchBounds {
    pub(crate) target_mz: f64,
    pub(crate) seed_rt: f64,
    pub(crate) rt_from: f64,
    pub(crate) rt_to: f64,
}

#[derive(Debug)]
pub(crate) struct GrowingCluster {
    items: Vec<TaggedFeature>,
    pub(crate) cached_median_mz: f64,
    pub(crate) cached_median_rt: f64,
}

impl GrowingCluster {
    pub(crate) fn new(item: TaggedFeature) -> Self {
        let mz = item.feature.mz;
        let rt = item.feature.rt;
        Self {
            items: vec![item],
            cached_median_mz: mz,
            cached_median_rt: rt,
        }
    }

    pub(crate) fn push(&mut self, item: TaggedFeature) {
        self.items.push(item);
        let mut mzs: Vec<f64> = self.items.iter().map(|t| t.feature.mz).collect();
        let mut rts: Vec<f64> = self.items.iter().map(|t| t.feature.rt).collect();
        mzs.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        rts.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        self.cached_median_mz = mzs[mzs.len() / 2];
        self.cached_median_rt = rts[rts.len() / 2];
    }

    fn into_items(self) -> Vec<TaggedFeature> {
        self.items
    }
}

impl FeatureClusterer {
    pub fn cluster(&self, mut tagged: Vec<TaggedFeature>) -> Vec<Cluster> {
        tagged.sort_unstable_by(|a, b| {
            a.feature
                .mz
                .partial_cmp(&b.feature.mz)
                .unwrap_or(Ordering::Equal)
        });
        let mz_groups = self.group_by_mz(tagged);
        mz_groups
            .into_iter()
            .flat_map(|group| self.subdivide_by_rt(group))
            .filter_map(|cluster| {
                let mut mzs: Vec<f64> = cluster.iter().map(|t| t.feature.mz).collect();
                let mut rts: Vec<f64> = cluster.iter().map(|t| t.feature.rt).collect();
                mzs.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
                rts.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
                let med_mz = mzs[mzs.len() / 2];
                let med_rt = rts[rts.len() / 2];
                let kept: Vec<TaggedFeature> = cluster
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

    fn group_by_mz(&self, items: Vec<TaggedFeature>) -> Vec<Vec<TaggedFeature>> {
        let mut groups: Vec<GrowingCluster> = Vec::new();
        for item in items {
            let belongs = groups.last().map_or(false, |g| {
                self.tolerance
                    .are_close_to_ref(item.feature.mz, g.cached_median_mz)
            });
            if belongs {
                groups.last_mut().unwrap().push(item);
            } else {
                groups.push(GrowingCluster::new(item));
            }
        }
        groups.into_iter().map(|g| g.into_items()).collect()
    }

    fn subdivide_by_rt(&self, group: Vec<TaggedFeature>) -> Vec<Vec<TaggedFeature>> {
        let mut rt_sorted = group;
        rt_sorted.sort_unstable_by(|a, b| {
            a.feature
                .rt
                .partial_cmp(&b.feature.rt)
                .unwrap_or(Ordering::Equal)
        });
        let mut clusters: Vec<GrowingCluster> = Vec::new();
        for item in rt_sorted {
            let belongs = clusters.last().map_or(false, |g| {
                (item.feature.rt - g.cached_median_rt).abs() <= self.rt_tolerance
            });
            if belongs {
                clusters.last_mut().unwrap().push(item);
            } else {
                clusters.push(GrowingCluster::new(item));
            }
        }
        clusters.into_iter().map(|g| g.into_items()).collect()
    }
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
pub fn get_features(
    directory_path: &str,
    time_window: FromTo,
    feature_config: FindFeaturesConfig,
    alignment_config: ConsensusAlignmentConfig,
    cores: usize,
) -> Result<Vec<ConsensusFeature>, AlignmentError> {
    let datasets = load_mzml_files(directory_path)?;
    let datasets = detect_features_per_sample(datasets, time_window, &feature_config, cores);

    let clusterer = FeatureClusterer {
        tolerance: alignment_config.tolerance.clone(),
        rt_tolerance: alignment_config.rt_tolerance,
    };

    let tagged: Vec<TaggedFeature> = datasets
        .iter()
        .enumerate()
        .flat_map(|(idx, (_, _, features))| {
            features.iter().map(move |f| TaggedFeature {
                sample_index: idx,
                feature: f.clone(),
            })
        })
        .collect();

    let clusters = clusterer.cluster(tagged);

    let results: Vec<ConsensusFeature> = clusters
        .into_iter()
        .filter_map(|cluster| build_consensus_feature(cluster, &datasets, &alignment_config))
        .collect();

    let points: Vec<(f64, f64, f64)> = results.iter().map(|f| (f.rt, f.mz, f.intensity)).collect();
    let survivors: HashSet<(u64, u64, u64)> = dedup_points(
        points.clone(),
        alignment_config.tolerance.mz_abs,
        alignment_config.tolerance.ppm,
        alignment_config.rt_tolerance,
    )
    .into_iter()
    .map(|(rt, mz, intensity)| (rt.to_bits(), mz.to_bits(), intensity.to_bits()))
    .collect();

    Ok(results
        .into_iter()
        .filter(|f| survivors.contains(&(f.rt.to_bits(), f.mz.to_bits(), f.intensity.to_bits())))
        .collect())
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
fn load_mzml_files(directory: &str) -> Result<Vec<(String, MzML)>, AlignmentError> {
    let entries: Vec<_> = fs::read_dir(directory)?
        .filter_map(|entry| {
            let path = entry.ok()?.path();
            matches!(path.extension()?.to_str()?, "mzML" | "b64").then_some(path)
        })
        .collect();

    if entries.is_empty() {
        return Err(AlignmentError::NoFiles);
    }

    entries
        .iter()
        .map(|path| {
            let file_name = path
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown")
                .to_string();
            let bytes = fs::read(path)?;
            let ext = path.extension().and_then(|s| s.to_str()).unwrap_or("");
            let mzml = match ext {
                "mzML" => parse_mzml(&bytes).map_err(|e| AlignmentError::Parse {
                    path: file_name.clone(),
                    source: e.to_string(),
                })?,
                "b64" => decode(&bytes).map_err(|e| AlignmentError::Parse {
                    path: file_name.clone(),
                    source: e.to_string(),
                })?,
                other => return Err(AlignmentError::UnsupportedFormat(other.to_string())),
            };
            Ok((file_name, mzml))
        })
        .collect()
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
fn detect_features_per_sample(
    samples: Vec<(String, MzML)>,
    time_window: FromTo,
    config: &FindFeaturesConfig,
    cores: usize,
) -> Vec<SampleDataset> {
    samples
        .into_iter()
        .enumerate()
        .map(|(idx, (name, mzml))| {
            let start = Instant::now();
            let features = find_features(
                &mzml,
                time_window,
                Some(FindFeaturesOptions {
                    scan_eic_options: Some(config.scan_eic_options),
                    eic_options: Some(config.eic_options),
                    find_peaks: Some(config.find_peaks.clone()),
                    mz_scan_grid: Some(config.mz_scan_grid.clone()),
                    scan_width_threshold: Some(config.scan_width_threshold),
                }),
                cores,
            )
            .unwrap_or_default();
            eprintln!(
                "[Sample {}] {} processed in {:.3}s",
                idx,
                name,
                start.elapsed().as_secs_f64()
            );
            (name, mzml, features)
        })
        .collect()
}

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
fn build_consensus_feature(
    cluster: Cluster,
    datasets: &[SampleDataset],
    config: &ConsensusAlignmentConfig,
) -> Option<ConsensusFeature> {
    let mut slots = assign_best_per_sample(cluster, datasets.len());
    let bounds = compute_search_bounds(&slots, config.rt_tolerance)?;
    fill_missing_slots(&mut slots, datasets, &bounds, config);
    let hits = collect_filled_slots(slots);
    require_minimum_frequency(hits, config.frequency)
        .map(|hits| aggregate_into_consensus(hits, &bounds))
}

pub(crate) fn assign_best_per_sample(cluster: Cluster, n_samples: usize) -> Vec<Option<Feature>> {
    let mut slots: Vec<Option<Feature>> = vec![None; n_samples];
    for tagged in cluster {
        let idx = tagged.sample_index;
        if slots[idx]
            .as_ref()
            .map_or(true, |f| tagged.feature.intensity > f.intensity)
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

#[cfg(not(all(target_arch = "wasm32", not(target_os = "wasi"))))]
fn fill_missing_slots(
    slots: &mut [Option<Feature>],
    datasets: &[SampleDataset],
    bounds: &SearchBounds,
    config: &ConsensusAlignmentConfig,
) {
    for (idx, slot) in slots.iter_mut().enumerate() {
        if slot.is_some() {
            continue;
        }
        if let Some(filled) = fill_missing_feature(
            &datasets[idx].1,
            bounds.target_mz,
            bounds.rt_from,
            bounds.rt_to,
            bounds.seed_rt,
            config.eic_options,
            config.peak_options.clone(),
        ) {
            if filled.intensity > 0.0 {
                *slot = Some(filled);
            }
        }
    }
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
        frequency: n,
    }
}

pub(crate) fn median(values: &mut [f64]) -> f64 {
    values.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
    let mid = values.len() / 2;
    if values.len() % 2 == 0 {
        (values[mid - 1] + values[mid]) / 2.0
    } else {
        values[mid]
    }
}

fn fill_missing_feature(
    mzml: &MzML,
    target_mz: f64,
    rt_start: f64,
    rt_end: f64,
    seed_rt: f64,
    eic_options: EicOptions,
    peak_options: Option<FindPeaksOptions>,
) -> Option<Feature> {
    let (time_points, ms1_scans) = collect_scans(
        mzml,
        FromTo {
            from: rt_start,
            to: rt_end,
        },
        TimeUnit::Minutes,
        1,
        false,
    );

    if ms1_scans.is_empty() {
        return None;
    }

    let intensities = compute_eic_for_mz(&ms1_scans, time_points.len(), target_mz, eic_options);
    let peak = get_peak(
        &DataXY {
            x: time_points,
            y: intensities,
        },
        Roi {
            rt: seed_rt,
            window: rt_end - rt_start,
        },
        peak_options,
    )?;

    if peak.intensity <= 0.0 {
        return None;
    }

    Some(Feature {
        mz: target_mz,
        rt: peak.rt,
        intensity: peak.intensity,
        from: peak.from,
        to: peak.to,
        np: peak.np,
        integral: peak.integral,
        noise: peak.noise,
    })
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
