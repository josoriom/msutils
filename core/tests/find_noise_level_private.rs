use msutils::utilities::{find_noise_level::find_noise_level, structs::DataXY};

use serde::Deserialize;
use std::{fs, path::Path};

mod helpers;
use helpers::approx_eq;

#[derive(Debug, Deserialize)]
struct Test {
    #[serde(alias = "metabolite", alias = "analyte")]
    name: String,
    data: DataXY,
    noise: f64,
}

fn load_tests() -> Vec<Test> {
    let path = Path::new(file!()).parent().unwrap().join("test.json");
    let s = fs::read_to_string(&path).unwrap_or_else(|e| panic!("cannot read {:?}: {}", path, e));
    serde_json::from_str(&s).expect("invalid JSON in tests/test.json")
}

fn find_test(name: &str) -> Test {
    let tests = load_tests();
    let test = tests
        .into_iter()
        .find(|json| json.name == name)
        .unwrap_or_else(|| panic!("entry not found in JSON: {}", name));
    return test;
}

#[test]
fn noise_matches_glutamic_acid() {
    let test = find_test("Glutamic acid");
    let noise = find_noise_level(&test.data.y) as f64;

    assert!(approx_eq(test.noise, noise, 0.1));
}

#[test]
fn noise_matches_proline_is() {
    let test = find_test("Proline IS");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_valine_is() {
    let test = find_test("Valine IS");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_aspartic_acid() {
    let test = find_test("Aspartic acid");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_leucine() {
    let test = find_test("Leucine");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_glutamine() {
    let test = find_test("Glutamine");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_alanine() {
    let test = find_test("Alanine");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_aminobutyric_acid() {
    let test = find_test("Aminobutyric acid");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}

#[test]
fn noise_matches_sarcosine() {
    let test = find_test("Sarcosine");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}

// ---- 10 Chrom tests ----
#[test]
fn noise_matches_chrom_1() {
    let test = find_test("Chrom 1");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_chrom_2() {
    let test = find_test("Chrom 2");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>>{noise:?}");
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_chrom_3() {
    let test = find_test("Chrom 3");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}

#[test]
fn noise_matches_chrom_4() {
    let test = find_test("Chrom 4");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_chrom_5() {
    let test = find_test("Chrom 5");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected: {}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_chrom_6() {
    let test = find_test("Chrom 6");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_chrom_7() {
    let test = find_test("Chrom 7");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_chrom_8() {
    let test = find_test("Chrom 8");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}
#[test]
fn noise_matches_chrom_9() {
    let test = find_test("Chrom 9");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}

#[test]
fn noise_matches_chrom_10() {
    let test = find_test("Chrom 10");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    assert!(approx_eq(test.noise, noise, 0.1));
}

#[test]
fn noise_matches_chrom_11() {
    let test = find_test("Chrom 11");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    print!("{noise}");
    assert!(approx_eq(test.noise, noise, 0.1));
}

#[test]
fn noise_matches_glutamic_acid_is() {
    let test = find_test("Glutamic acid - IS");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    print!("{noise}");
    assert!(approx_eq(test.noise, noise, 0.1));
}

#[test]
fn noise_matches_beta_alanine() {
    let test = find_test("beta-alanine");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    print!("{noise}");
    assert!(approx_eq(test.noise, noise, 0.1));
}

#[test]
fn noise_matches_taurine() {
    let test = find_test("Taurine");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    print!("{noise}");
    assert!(approx_eq(test.noise, noise, 0.1));
}

#[test]
#[ignore = "I'll check later"]
fn noise_matches_chrom_12() {
    let test = find_test("Chrom 12");
    let noise = find_noise_level(&test.data.y) as f64;
    println!("---::>> got: {}, expected{}", noise, test.noise);
    print!("{noise}");
    assert!(approx_eq(test.noise, noise, 0.1));
}
