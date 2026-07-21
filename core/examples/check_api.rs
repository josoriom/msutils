use std::path::Path;

use ionic::ion::{IonReader, ReadOptions};
use quantion::utilities::{
    calculate_eic::{EicOptions, EicReader, calculate_eic, get_scan_times},
    find_features::{FindFeaturesOptions, find_features},
    find_noise_level::find_noise_level,
    find_peaks::{FindPeaksOptions, PeakFilter, find_peaks},
    get_features::{AlignmentOptions, get_features},
    get_peak::get_peak,
    get_peaks_from_chrom::get_peaks_from_chrom,
    get_peaks_from_eic::get_peaks_from_eic,
    structs::{DataXY, FromTo, Peak, Roi},
};

const FROM_RT: f64 = 0.0;
const TO_RT: f64 = 5.0;
const EIC_PPM: f64 = 20.0;
const EIC_MZ: f64 = 0.005;
const GRID_STEP: f64 = 0.005;
const GRID_START: f64 = 89.0;
const GRID_END: f64 = 91.0;
const ROI_RANGE: f64 = 0.05;
const CORES: usize = 1;
const MS1_LEVEL: u8 = 1;

const FEATURE_MIN_INTENSITY: f64 = 500.0;
const FEATURE_MIN_PEAK_WIDTH_POINTS: usize = 5;

const SELF_CHECK_VALUES: [f64; 3] = [1.0, 2.0, 3.0];

struct Channel {
    name: &'static str,
    mass: f64,
    min_intensity: f64,
    min_peak_width_points: usize,
    targets: [(&'static str, f64); 3],
}

const CHANNELS: [Channel; 2] = [
    Channel {
        name: "mz_2",
        mass: 90.05550,
        min_intensity: 500.0,
        min_peak_width_points: 5,
        targets: [("met_1", 2.885), ("met_2", 2.552), ("met_3", 2.465)],
    },
    Channel {
        name: "mz_1",
        mass: 89.04768,
        min_intensity: 500.0,
        min_peak_width_points: 5,
        targets: [("met_1", 2.385), ("met_2", 2.18), ("met_3", 2.08)],
    },
];

fn number(value: f64) -> String {
    if value.is_nan() {
        return "nan".to_string();
    }
    if value.is_infinite() {
        return if value > 0.0 { "inf" } else { "-inf" }.to_string();
    }
    format!("{value}")
}

fn digest(values: &[f64]) -> String {
    let mut hash: u32 = 0x811c_9dc5;
    for value in values {
        for byte in value.to_bits().to_le_bytes() {
            hash ^= byte as u32;
            hash = hash.wrapping_mul(0x0100_0193);
        }
    }
    format!("{hash:08x}")
}

fn eic_options() -> EicOptions {
    EicOptions {
        ppm_tolerance: EIC_PPM,
        mz_tolerance: EIC_MZ,
        ..Default::default()
    }
}
fn peak_options(channel: &Channel) -> FindPeaksOptions {
    FindPeaksOptions {
        filter: Some(PeakFilter {
            min_intensity: Some(channel.min_intensity),
            min_peak_width_points: Some(channel.min_peak_width_points),
            auto_noise: Some(true),
            auto_baseline: Some(true),
            ..Default::default()
        }),
        ..Default::default()
    }
}

fn chrom_peak_options(channel: &Channel) -> FindPeaksOptions {
    FindPeaksOptions {
        filter: Some(PeakFilter {
            min_intensity: Some(channel.min_intensity),
            min_peak_width_points: Some(channel.min_peak_width_points),
            auto_noise: Some(false),
            auto_baseline: Some(false),
            ..Default::default()
        }),
        ..Default::default()
    }
}

fn feature_peak_options() -> FindPeaksOptions {
    FindPeaksOptions {
        filter: Some(PeakFilter {
            min_intensity: Some(FEATURE_MIN_INTENSITY),
            min_peak_width_points: Some(FEATURE_MIN_PEAK_WIDTH_POINTS),
            auto_noise: Some(true),
            auto_baseline: Some(true),
            ..Default::default()
        }),
        ..Default::default()
    }
}

fn feature_options() -> FindFeaturesOptions {
    let mut options = FindFeaturesOptions::default();
    options.final_eic_options.ppm_tolerance = EIC_PPM;
    options.final_eic_options.mz_tolerance = EIC_MZ;
    options.mz_scan_grid.min_mz = GRID_START;
    options.mz_scan_grid.max_mz = GRID_END;
    options.mz_scan_grid.step = GRID_STEP;
    options.peak_options = feature_peak_options();
    options
}

fn alignment_options() -> AlignmentOptions {
    let mut options = AlignmentOptions::default();
    options.eic_options = eic_options();
    options.peak_options = Some(feature_peak_options());
    options.min_samples = 1;
    options
}

fn time_range() -> FromTo {
    FromTo {
        from: FROM_RT,
        to: TO_RT,
    }
}

fn open(fixture: &Path) -> IonReader {
    IonReader::open_file(fixture, ReadOptions::default()).expect("open fixture")
}

fn count_scans(fixture: &Path) -> usize {
    let mut ion = open(fixture);
    let mut reader = EicReader::Ion(&mut ion);
    get_scan_times(&mut reader, FROM_RT, TO_RT, MS1_LEVEL).len()
}

fn extract_eic(fixture: &Path, channel: &Channel) -> DataXY {
    let mut ion = open(fixture);
    let mut reader = EicReader::Ion(&mut ion);
    calculate_eic(&mut reader, channel.mass, time_range(), eic_options()).expect("calculate_eic")
}

fn peak_lines(prefix: &str, peak: &Peak, out: &mut Vec<String>) {
    shared_peak_lines(prefix, peak, out);
    out.push(format!("{prefix}.n_points = {}", peak.n_points));
}

fn shared_peak_lines(prefix: &str, peak: &Peak, out: &mut Vec<String>) {
    out.push(format!("{prefix}.rt = {}", number(peak.rt)));
    out.push(format!("{prefix}.from = {}", number(peak.from)));
    out.push(format!("{prefix}.to = {}", number(peak.to)));
    out.push(format!("{prefix}.integral = {}", number(peak.integral)));
    out.push(format!("{prefix}.intensity = {}", number(peak.intensity)));
}

fn peaks_from_eic(fixture: &Path, channel: &Channel, prefix: &str, out: &mut Vec<String>) {
    let mut ion = open(fixture);
    let mut reader = EicReader::Ion(&mut ion);
    let rois: Vec<Roi> = channel
        .targets
        .iter()
        .map(|(id, rt)| Roi::eic(*id, *rt, channel.mass, ROI_RANGE))
        .collect();
    let rows = get_peaks_from_eic(
        &mut reader,
        time_range(),
        &rois,
        Some(peak_options(channel)),
        CORES,
    )
    .expect("get_peaks_from_eic");
    for (id, _, _, peak) in &rows {
        shared_peak_lines(&format!("{prefix}.{id}"), peak, out);
    }
}

fn chromatogram_index(channel: &Channel) -> usize {
    CHANNELS
        .iter()
        .filter(|other| other.mass < channel.mass)
        .count()
}

fn peaks_from_chrom(fixture: &Path, channel: &Channel, prefix: &str, out: &mut Vec<String>) {
    let mut ion = open(fixture);
    let mzml = ion.to_mzml().expect("read chromatograms");
    let index = chromatogram_index(channel);
    let rois: Vec<Roi> = channel
        .targets
        .iter()
        .map(|(id, rt)| Roi::chrom(*id, index, *rt, ROI_RANGE))
        .collect();
    let rows = get_peaks_from_chrom(&mzml, &rois, Some(chrom_peak_options(channel)), CORES)
        .unwrap_or_default();
    for (roi, row) in channel.targets.iter().zip(&rows) {
        let name = format!("{prefix}.{}", roi.0);
        out.push(format!("{name}.rt = {}", number(row.peak_rt)));
        out.push(format!("{name}.from = {}", number(row.from_rt)));
        out.push(format!("{name}.to = {}", number(row.to_rt)));
        out.push(format!("{name}.intensity = {}", number(row.intensity)));
        out.push(format!("{name}.area = {}", number(row.area)));
    }
}

fn features(fixture: &Path, out: &mut Vec<String>) {
    let mut ion = open(fixture);
    let mut reader = EicReader::Ion(&mut ion);
    let found = find_features(&mut reader, time_range(), Some(feature_options()), CORES)
        .expect("find_features");
    out.push(format!("find_features.count = {}", found.len()));
    for (index, f) in found.iter().enumerate() {
        let name = format!("find_features[{index}]");
        out.push(format!("{name}.mz = {}", number(f.mz)));
        out.push(format!("{name}.rt = {}", number(f.rt)));
        out.push(format!("{name}.from = {}", number(f.from)));
        out.push(format!("{name}.to = {}", number(f.to)));
        out.push(format!("{name}.integral = {}", number(f.integral)));
        out.push(format!("{name}.intensity = {}", number(f.intensity)));
        out.push(format!("{name}.n_points = {}", f.n_points));
    }
}

fn consensus_features(directory: &str, out: &mut Vec<String>) {
    let found = get_features(
        directory,
        time_range(),
        feature_options(),
        alignment_options(),
        CORES,
    )
    .expect("get_features");
    out.push(format!("get_features.count = {}", found.len()));
    for (index, f) in found.iter().enumerate() {
        let name = format!("get_features[{index}]");
        out.push(format!("{name}.mz = {}", number(f.mz)));
        out.push(format!("{name}.rt = {}", number(f.rt)));
        out.push(format!("{name}.from = {}", number(f.from)));
        out.push(format!("{name}.to = {}", number(f.to)));
        out.push(format!("{name}.integral = {}", number(f.integral)));
        out.push(format!("{name}.intensity = {}", number(f.intensity)));
    }
}

fn report_channel(fixture: &Path, channel: &Channel, out: &mut Vec<String>) {
    let name = channel.name;
    let eic = extract_eic(fixture, channel);
    let noise = find_noise_level(eic.y.as_slice());
    let data = DataXY {
        x: eic.x.clone(),
        y: eic.y.clone(),
    };

    out.push(format!("{name}.calculate_eic.points = {}", eic.x.len()));
    out.push(format!("{name}.calculate_eic.x.fnv1a = {}", digest(&eic.x)));
    out.push(format!("{name}.calculate_eic.y.fnv1a = {}", digest(&eic.y)));
    out.push(format!(
        "{name}.find_noise_level.intensity = {}",
        number(noise.intensity)
    ));
    out.push(format!("{name}.find_noise_level.width = {}", noise.width));

    let found = find_peaks(&data, Some(peak_options(channel)));
    out.push(format!("{name}.find_peaks.count = {}", found.len()));
    for (index, peak) in found.iter().enumerate() {
        peak_lines(&format!("{name}.find_peaks[{index}]"), peak, out);
    }

    let target = channel.targets[2].1;
    let single = get_peak(
        &data,
        &Roi::peak(target, ROI_RANGE),
        Some(peak_options(channel)),
    );
    peak_lines(&format!("{name}.get_peak"), &single, out);

    peaks_from_eic(fixture, channel, &format!("{name}.get_peaks_from_eic"), out);
    peaks_from_chrom(
        fixture,
        channel,
        &format!("{name}.get_peaks_from_chrom"),
        out,
    );
}

fn main() {
    let directory = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "core/tests/fixtures/api".to_string());
    let fixture = Path::new(&directory).join("api.ion");

    let mut out = Vec::new();
    out.push(format!("fnv1a.self_check = {}", digest(&SELF_CHECK_VALUES)));
    out.push(format!("parse_ion.scans = {}", count_scans(&fixture)));
    for channel in CHANNELS.iter() {
        report_channel(&fixture, channel, &mut out);
    }
    features(&fixture, &mut out);
    consensus_features(&directory, &mut out);

    println!("{}", out.join("\n"));
}
