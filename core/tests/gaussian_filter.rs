use msutils::utilities::find_peaks::{ArtifactFilter, FindPeaksOptions, PeakFilter};
use msutils::utilities::get_peak::get_peak;
use msutils::utilities::structs::{DataXY, Roi};

use serde::Deserialize;
use std::{collections::HashMap, fs, path::Path};

mod helpers;

#[derive(Debug, Deserialize)]
struct Feature {
    mz: f64,
    rt: f64,
    x: Vec<f64>,
    y: Vec<f64>,
    #[serde(rename = "isPeak")]
    is_peak: bool,
}

fn load_features() -> HashMap<String, Feature> {
    let path = Path::new(file!())
        .parent()
        .unwrap()
        .join("gaussian_filter.json");
    let s = fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {:?}: {}", path, e));
    serde_json::from_str(&s).expect("invalid JSON in gaussian_filter.json")
}

fn default_options() -> FindPeaksOptions {
    FindPeaksOptions {
        filter: Some(PeakFilter {
            auto_noise: Some(true),
            auto_baseline: Some(true),
            min_peak_width_points: Some(5),
            min_snr: Some(2.0),
            ..Default::default()
        }),
        artifact_filter: Some(ArtifactFilter {
            min_gaussian_r2: 0.30,
            min_apex_to_boundary_ratio: 2.0,
        }),
        ..Default::default()
    }
}

fn peaks_near_rt(data: &DataXY, rt: f64, opts: FindPeaksOptions) -> bool {
    let roi = Roi::new(rt, 0.05);
    let peak = get_peak(data, &roi, Some(opts));
    if let Some(p) = peak {
        helpers::dump_peaks(&[p]);
        (p.rt - rt).abs() < 0.05 && p.intensity > 0.0
    } else {
        false
    }
}

#[test]
fn feature_1_mz42_rt0617_is_not_a_peak() {
    let features = load_features();
    let f = &features["feature_1"];
    let data = DataXY {
        x: f.x.clone(),
        y: f.y.clone(),
    };
    assert!(!f.is_peak);
    assert!(
        !peaks_near_rt(&data, f.rt, default_options()),
        "mz={:.4} rt={:.4} should be filtered as artifact",
        f.mz,
        f.rt
    );
}

#[test]
fn feature_2_mz43_rt0637_is_not_a_peak() {
    let features = load_features();
    let f = &features["feature_2"];
    let data = DataXY {
        x: f.x.clone(),
        y: f.y.clone(),
    };
    assert!(!f.is_peak);
    assert!(
        !peaks_near_rt(&data, f.rt, default_options()),
        "mz={:.4} rt={:.4} should be filtered as artifact",
        f.mz,
        f.rt
    );
}

#[test]
fn feature_3_mz61_rt060_is_not_a_peak() {
    let features = load_features();
    let f = &features["feature_3"];
    let data = DataXY {
        x: f.x.clone(),
        y: f.y.clone(),
    };
    assert!(!f.is_peak);
    assert!(
        !peaks_near_rt(&data, f.rt, default_options()),
        "mz={:.4} rt={:.4} should be filtered as artifact",
        f.mz,
        f.rt
    );
}

#[test]
fn feature_4_mz281_rt070_is_not_a_peak() {
    let features = load_features();
    let f = &features["feature_4"];
    let data = DataXY {
        x: f.x.clone(),
        y: f.y.clone(),
    };
    assert!(!f.is_peak);
    assert!(
        !peaks_near_rt(&data, f.rt, default_options()),
        "mz={:.4} rt={:.4} should be filtered as artifact",
        f.mz,
        f.rt
    );
}

#[test]
fn feature_5_mz203_rt0716_is_not_a_peak() {
    let features = load_features();
    let f = &features["feature_5"];
    let data = DataXY {
        x: f.x.clone(),
        y: f.y.clone(),
    };
    assert!(!f.is_peak);
    assert!(
        !peaks_near_rt(&data, f.rt, default_options()),
        "mz={:.4} rt={:.4} should be filtered as artifact",
        f.mz,
        f.rt
    );
}

#[test]
fn feature_6_mz90_rt061_is_not_a_peak() {
    let features = load_features();
    let f = &features["feature_6"];
    let data = DataXY {
        x: f.x.clone(),
        y: f.y.clone(),
    };
    assert!(!f.is_peak);
    assert!(
        !peaks_near_rt(&data, f.rt, default_options()),
        "mz={:.4} rt={:.4} should be filtered as artifact",
        f.mz,
        f.rt
    );
}

#[ignore = "I'll check later"]
#[test]
fn feature_7_mz84_rt0567_is_a_peak() {
    let features = load_features();
    let f = &features["feature_7"];
    let data = DataXY {
        x: f.x.clone(),
        y: f.y.clone(),
    };
    assert!(f.is_peak);
    assert!(
        peaks_near_rt(&data, f.rt, default_options()),
        "mz={:.4} rt={:.4} is a real peak and must not be filtered",
        f.mz,
        f.rt
    );
}

#[test]
fn feature_8_mz279_rt0569_is_a_peak() {
    let features = load_features();
    let f = &features["feature_8"];
    let data = DataXY {
        x: f.x.clone(),
        y: f.y.clone(),
    };
    assert!(f.is_peak);
    assert!(
        peaks_near_rt(&data, f.rt, default_options()),
        "mz={:.4} rt={:.4} is a real peak and must not be filtered",
        f.mz,
        f.rt
    );
}

#[test]
fn feature_9_mz808_rt0661_is_a_peak() {
    let features = load_features();
    let f = &features["feature_9"];
    let data = DataXY {
        x: f.x.clone(),
        y: f.y.clone(),
    };
    assert!(f.is_peak);
    assert!(
        peaks_near_rt(&data, f.rt, default_options()),
        "mz={:.4} rt={:.4} is a real peak and must not be filtered",
        f.mz,
        f.rt
    );
}

#[test]
fn feature_10() {
    let features = load_features();
    let f = &features["feature_10"];
    let data = DataXY {
        x: f.x.clone(),
        y: f.y.clone(),
    };
    let should_pass = f.is_peak;
    let found = peaks_near_rt(&data, f.rt, default_options());
    if should_pass {
        assert!(
            found,
            "mz={:.4} rt={:.4} is a real peak and must not be filtered",
            f.mz, f.rt
        );
    } else {
        assert!(
            !found,
            "mz={:.4} rt={:.4} should be filtered as artifact",
            f.mz, f.rt
        );
    }
}

#[test]
fn feature_11() {
    let features = load_features();
    let f = &features["feature_11"];
    let data = DataXY {
        x: f.x.clone(),
        y: f.y.clone(),
    };
    let should_pass = f.is_peak;
    let found = peaks_near_rt(&data, f.rt, default_options());
    if should_pass {
        assert!(
            found,
            "mz={:.4} rt={:.4} is a real peak and must not be filtered",
            f.mz, f.rt
        );
    } else {
        assert!(
            !found,
            "mz={:.4} rt={:.4} should be filtered as artifact",
            f.mz, f.rt
        );
    }
}

#[test]
fn feature_12() {
    let features = load_features();
    let f = &features["feature_12"];
    let data = DataXY {
        x: f.x.clone(),
        y: f.y.clone(),
    };
    let should_pass = f.is_peak;
    let found = peaks_near_rt(&data, f.rt, default_options());
    if should_pass {
        assert!(
            found,
            "mz={:.4} rt={:.4} is a real peak and must not be filtered",
            f.mz, f.rt
        );
    } else {
        assert!(
            !found,
            "mz={:.4} rt={:.4} should be filtered as artifact",
            f.mz, f.rt
        );
    }
}

#[test]
fn feature_13() {
    let features = load_features();
    let f = &features["feature_13"];
    let data = DataXY {
        x: f.x.clone(),
        y: f.y.clone(),
    };
    let should_pass = f.is_peak;
    let found = peaks_near_rt(&data, f.rt, default_options());
    if should_pass {
        assert!(
            found,
            "mz={:.4} rt={:.4} is a real peak and must not be filtered",
            f.mz, f.rt
        );
    } else {
        assert!(
            !found,
            "mz={:.4} rt={:.4} should be filtered as artifact",
            f.mz, f.rt
        );
    }
}

#[ignore = "I'll check later"]
#[test]
fn feature_14() {
    let features = load_features();
    let f = &features["feature_14"];
    let data = DataXY {
        x: f.x.clone(),
        y: f.y.clone(),
    };
    let should_pass = f.is_peak;
    let found = peaks_near_rt(&data, f.rt, default_options());
    if should_pass {
        assert!(
            found,
            "mz={:.4} rt={:.4} is a real peak and must not be filtered",
            f.mz, f.rt
        );
    } else {
        assert!(
            !found,
            "mz={:.4} rt={:.4} should be filtered as artifact",
            f.mz, f.rt
        );
    }
}

#[test]
fn feature_15() {
    let features = load_features();
    let f = &features["feature_15"];
    let data = DataXY {
        x: f.x.clone(),
        y: f.y.clone(),
    };
    let should_pass = f.is_peak;
    let found = peaks_near_rt(&data, f.rt, default_options());
    if should_pass {
        assert!(
            found,
            "mz={:.4} rt={:.4} is a real peak and must not be filtered",
            f.mz, f.rt
        );
    } else {
        assert!(
            !found,
            "mz={:.4} rt={:.4} should be filtered as artifact",
            f.mz, f.rt
        );
    }
}

#[test]
fn feature_16() {
    let features = load_features();
    let f = &features["feature_16"];
    let data = DataXY {
        x: f.x.clone(),
        y: f.y.clone(),
    };
    let should_pass = f.is_peak;
    let found = peaks_near_rt(&data, f.rt, default_options());
    println!("feature_16: isPeak={}, found={}", should_pass, found);
    if should_pass {
        assert!(
            found,
            "mz={:.4} rt={:.4} is a real peak and must not be filtered",
            f.mz, f.rt
        );
    } else {
        assert!(
            !found,
            "mz={:.4} rt={:.4} should be filtered as artifact",
            f.mz, f.rt
        );
    }
}

#[test]
fn feature_17() {
    let features = load_features();
    let f = &features["feature_17"];
    let data = DataXY {
        x: f.x.clone(),
        y: f.y.clone(),
    };
    let should_pass = f.is_peak;
    let found = peaks_near_rt(&data, f.rt, default_options());
    println!("feature_17: isPeak={}, found={}", should_pass, found);
    if should_pass {
        assert!(
            found,
            "mz={:.4} rt={:.4} is a real peak and must not be filtered",
            f.mz, f.rt
        );
    } else {
        assert!(
            !found,
            "mz={:.4} rt={:.4} should be filtered as artifact",
            f.mz, f.rt
        );
    }
}

#[ignore = "I'll check later"]
#[test]
fn feature_18() {
    let features = load_features();
    let f = &features["feature_18"];
    let data = DataXY {
        x: f.x.clone(),
        y: f.y.clone(),
    };
    let should_pass = f.is_peak;
    let found = peaks_near_rt(&data, f.rt, default_options());
    println!("feature_18: isPeak={}, found={}", should_pass, found);
    if should_pass {
        assert!(
            found,
            "mz={:.4} rt={:.4} is a real peak and must not be filtered",
            f.mz, f.rt
        );
    } else {
        assert!(
            !found,
            "mz={:.4} rt={:.4} should be filtered as artifact",
            f.mz, f.rt
        );
    }
}
