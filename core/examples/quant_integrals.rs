use std::collections::HashSet;
use std::env;
use std::fs;
use std::io::Write;

use ionic::ion::{IonReader, ReadOptions};
use msutils::utilities::calculate_eic::{Eic, EicOptions, EicReader, calculate_eic};
use msutils::utilities::find_peaks::{FindPeaksOptions, PeakFilter};
use msutils::utilities::get_peak::get_peak;
use msutils::utilities::structs::{DataXY, FromTo, Roi};
use serde_json::{Value, json};

const HALF_WIDTH: f64 = 0.5;
const EIC_WINDOW: f64 = 1.0;

fn options() -> FindPeaksOptions {
    FindPeaksOptions {
        filter: Some(PeakFilter {
            auto_noise: Some(true),
            auto_baseline: Some(true),
            min_peak_width_points: Some(5),
            min_intensity: Some(1000.0),
            min_snr: Some(3.0),
            ..Default::default()
        }),
        ..Default::default()
    }
}

fn nearest_index(x: &[f64], value: f64) -> Option<usize> {
    if x.is_empty() || !value.is_finite() {
        return None;
    }
    let mut best = 0usize;
    let mut best_distance = f64::INFINITY;
    for (index, point) in x.iter().enumerate() {
        let distance = (point - value).abs();
        if distance < best_distance {
            best_distance = distance;
            best = index;
        }
    }
    Some(best)
}

fn column(header: &[&str], name: &str) -> usize {
    header
        .iter()
        .position(|field| *field == name)
        .unwrap_or_else(|| panic!("missing column {name}"))
}

fn read_reference(path: &str) -> (Vec<(String, f64, f64)>, Vec<String>) {
    let text = fs::read_to_string(path).expect("read reference tsv");
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().expect("reference header").split('\t').collect();
    let sample_column = column(&header, "sample");
    let id_column = column(&header, "id");
    let mz_column = column(&header, "mz");
    let rt_column = column(&header, "rt");

    let mut features = Vec::new();
    let mut seen_features = HashSet::new();
    let mut samples = Vec::new();
    let mut seen_samples = HashSet::new();

    for line in lines {
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let id = fields[id_column].to_string();
        if seen_features.insert(id.clone()) {
            let mz = fields[mz_column].parse::<f64>().unwrap_or(f64::NAN);
            let rt = fields[rt_column].parse::<f64>().unwrap_or(f64::NAN);
            features.push((id, mz, rt));
        }
        let sample = fields[sample_column].to_string();
        if seen_samples.insert(sample.clone()) {
            samples.push(sample);
        }
    }
    (features, samples)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let reference = &args[1];
    let ion_dir = &args[2];
    let output = &args[3];
    let limit = args.get(4).and_then(|value| value.parse::<usize>().ok());

    let (mut features, samples) = read_reference(reference);
    if let Some(count) = limit {
        features.truncate(count);
    }

    let mut eics_by_feature: Vec<Vec<Value>> = vec![Vec::new(); features.len()];
    let mut rows = vec![String::from("id\tsample\tfrom_rt\tto_rt\tinteg\tn_points")];

    for sample in &samples {
        let path = format!("{ion_dir}/{sample}.ion");
        let bytes = fs::read(&path).expect("read ion");
        let mut ion = IonReader::open(&bytes, ReadOptions::default()).expect("open ion");
        let mut reader = EicReader::Ion(&mut ion);
        eprintln!("processing {sample}...");

        for (index, (id, mz, rt)) in features.iter().enumerate() {
            let eic = match calculate_eic(
                &mut reader,
                *mz,
                FromTo {
                    from: rt - EIC_WINDOW,
                    to: rt + EIC_WINDOW,
                },
                EicOptions::default(),
            ) {
                Ok(value) => value,
                Err(_) => Eic {
                    x: Vec::new(),
                    y: Vec::new(),
                },
            };

            let data = DataXY {
                x: eic.x,
                y: eic.y,
            };
            let peak = get_peak(&data, &Roi::new(*rt, HALF_WIDTH), Some(options()));

            let eic_entry = match &peak {
                Some(found) => {
                    let from_index = nearest_index(&data.x, found.from);
                    let to_index = nearest_index(&data.x, found.to);
                    rows.push(format!(
                        "{id}\t{sample}\t{}\t{}\t{}\t{}",
                        found.from, found.to, found.integral, found.n_points
                    ));
                    json!({
                        "sample": sample,
                        "x": data.x,
                        "y": data.y,
                        "from": { "rt": found.from, "index": from_index },
                        "to": { "rt": found.to, "index": to_index },
                        "integral": found.integral,
                    })
                }
                None => {
                    rows.push(format!("{id}\t{sample}\tNA\tNA\t0\t0"));
                    json!({
                        "sample": sample,
                        "x": data.x,
                        "y": data.y,
                        "from": { "rt": Value::Null, "index": Value::Null },
                        "to": { "rt": Value::Null, "index": Value::Null },
                        "integral": Value::Null,
                    })
                }
            };
            eics_by_feature[index].push(eic_entry);
        }
    }

    let mut features_json = Vec::with_capacity(features.len());
    for (index, (id, mz, rt)) in features.iter().enumerate() {
        features_json.push(json!({
            "id": id,
            "mz": mz,
            "rt": rt,
            "eics": std::mem::take(&mut eics_by_feature[index]),
        }));
    }

    let serialized = serde_json::to_string(&Value::Array(features_json)).expect("serialize output");
    let mut file = fs::File::create(output).expect("create output");
    file.write_all(serialized.as_bytes()).expect("write output");

    print!("{}", rows.join("\n"));
    println!();
    eprintln!(
        "wrote {output} ({} features x {} samples)",
        features.len(),
        samples.len()
    );
}
