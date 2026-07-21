use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
    sync::OnceLock,
};

use quantion::utilities::{find_noise_level::find_noise_level, structs::DataXY};

mod helpers;
use helpers::{Eic, approx_eq, load_chromatograms_from};

struct Test {
    data: DataXY,
    noise: f64,
}

fn ion_path() -> PathBuf {
    if let Ok(path) = std::env::var("TEST_ION") {
        return PathBuf::from(path);
    }
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("fixtures")
        .join("test.ion")
}

fn all_cases() -> &'static BTreeMap<String, Eic> {
    static CACHE: OnceLock<BTreeMap<String, Eic>> = OnceLock::new();
    CACHE.get_or_init(|| load_chromatograms_from(&ion_path()))
}

fn find_test(name: &str) -> Test {
    let eic = all_cases()
        .get(name)
        .unwrap_or_else(|| panic!("noise case not found in test.ion: {name}"));
    let noise = eic
        .params
        .get("exp.noise")
        .unwrap_or_else(|| panic!("case {name} has no exp.noise"))
        .parse()
        .expect("numeric exp.noise");
    Test {
        data: DataXY {
            x: eic.time.clone(),
            y: eic.intensity.clone(),
        },
        noise,
    }
}

#[test]
fn noise_matches_glutamic_acid() {
    let test = find_test("Glutamic acid");
    let noise = find_noise_level(&test.data.y).intensity as f64;

    assert!(approx_eq(test.noise, noise, 0.1));
}

#[test]
fn noise_matches_proline_is() {
    let test = find_test("Proline IS");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_valine_is() {
    let test = find_test("Valine IS");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_aspartic_acid() {
    let test = find_test("Aspartic acid");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_leucine() {
    let test = find_test("Leucine");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_glutamine() {
    let test = find_test("Glutamine");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_alanine() {
    let test = find_test("Alanine");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_aminobutyric_acid() {
    let test = find_test("Aminobutyric acid");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}

#[test]
fn noise_matches_sarcosine() {
    let test = find_test("Sarcosine");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}

// ---- 10 Chrom tests ----
#[test]
fn noise_matches_chrom_1() {
    let test = find_test("Chrom 1");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_chrom_2() {
    let test = find_test("Chrom 2");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>>{noise:?}");
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_chrom_3() {
    let test = find_test("Chrom 3");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}

#[test]
fn noise_matches_chrom_4() {
    let test = find_test("Chrom 4");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_chrom_5() {
    let test = find_test("Chrom 5");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected: {}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_chrom_6() {
    let test = find_test("Chrom 6");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_chrom_7() {
    let test = find_test("Chrom 7");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_chrom_8() {
    let test = find_test("Chrom 8");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_chrom_9() {
    let test = find_test("Chrom 9");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}

#[test]
fn noise_matches_chrom_10() {
    let test = find_test("Chrom 10");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}

#[test]
#[ignore = "I'll check later"]
fn noise_matches_chrom_11() {
    let test = find_test("Chrom 11");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    print!("{noise}");
    assert!(approx_eq(test.noise, noise, 0.1));
}

#[test]
fn noise_matches_glutamic_acid_is() {
    let test = find_test("Glutamic acid - IS");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    print!("{noise}");
    assert!(approx_eq(test.noise, noise, 0.1));
}

#[test]
fn noise_matches_beta_alanine() {
    let test = find_test("beta-alanine");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    print!("{noise}");
    assert!(approx_eq(test.noise, noise, 0.1));
}

#[test]
fn noise_matches_taurine() {
    let test = find_test("Taurine");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    print!("{noise}");
    assert!(approx_eq(test.noise, noise, 0.1));
}

#[test]
#[ignore = "I'll check later"]
fn noise_matches_chrom_12() {
    let test = find_test("Chrom 12");
    let noise = find_noise_level(&test.data.y).intensity as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    print!("{noise}");
    assert!(approx_eq(test.noise, noise, 0.1));
}
