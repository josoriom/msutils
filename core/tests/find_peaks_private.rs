use msutils::utilities::calculate_baseline::BaselineOptions;
use msutils::utilities::find_peaks::{FilterPeaksOptions, FindPeaksOptions, find_peaks};
use msutils::utilities::structs::DataXY;

use serde::Deserialize;
use std::{fs, path::Path};

mod helpers;
use helpers::approx_eq;

use crate::helpers::dump_peaks;

#[derive(Debug, Deserialize)]
struct Test {
    #[serde(alias = "metabolite", alias = "analyte")]
    name: String,
    data: DataXY,
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
fn find_3_peaks_aminobutyric_acid_1() {
    let test = find_test("Aminobutyric acid");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(10),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);

    assert_eq!(peaks.len(), 2);

    // RT
    assert!(approx_eq(peaks[0].rt, 2.55, 0.06));
    assert!(approx_eq(peaks[1].rt, 2.89, 0.06));
    // Intensity
    assert!(approx_eq(peaks[0].intensity, 1204.0, 1.0));
    assert!(approx_eq(peaks[1].intensity, 30580.0, 1.0));
    // Integral
    // assert!(approx_eq(peaks[0].integral, 28.345, 5.0));
    // assert!(approx_eq(peaks[1].integral, 721.216, 10.0));
}

#[test]
fn find_3_peaks_alanine_acid_1() {
    let test = find_test("Alanine");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(5),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);

    assert_eq!(peaks.len(), 3);

    // RT
    assert!(approx_eq(peaks[0].rt, 2.09, 0.06));
    assert!(approx_eq(peaks[1].rt, 2.19, 0.06));
    assert!(approx_eq(peaks[2].rt, 2.39, 0.06));
    // Intensity
    assert!(approx_eq(peaks[0].intensity, 1086.0, 1.0));
    assert!(approx_eq(peaks[1].intensity, 2974.0, 1.0));
    assert!(approx_eq(peaks[2].intensity, 157154.0, 1.0));
    // Integral
    // assert!(approx_eq(peaks[0].integral, 27.32, 10.0));
    // assert!(approx_eq(peaks[1].integral, 73.84, 10.0));
    // assert!(approx_eq(peaks[2].integral, 4238.05, 100.0));
}

#[test]
fn find_1_peak_glutaime_with_tail() {
    let test = find_test("Glutamine");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_baseline: Some(true),
                auto_noise: Some(true),
                width_threshold: Some(8),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);

    assert_eq!(peaks.len(), 1);

    // RT
    assert!(approx_eq(peaks[0].rt, 1.91, 0.05));
    // Intensity
    assert!(approx_eq(peaks[0].intensity, 221458.00, 100.0));
    // Integral
    // assert!(approx_eq(peaks[0].integral, 6219.72, 100.0));
}

#[test]
fn find_2_peaks_leucine() {
    let test = find_test("Leucine");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(8),
                intensity_threshold: Some(400.0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);

    assert_eq!(peaks.len(), 2);

    // RT
    assert!(approx_eq(peaks[0].rt, 4.07, 0.05));
    assert!(approx_eq(peaks[1].rt, 4.15, 0.05));
    // Intensity
    assert!(approx_eq(peaks[0].intensity, 43570.0, 100.0));
    assert!(approx_eq(peaks[1].intensity, 90706.00, 100.0));
    // Integral
    // assert!(approx_eq(peaks[0].integral, 955.969, 10.0));
    // assert!(approx_eq(peaks[1].integral, 2359.649, 20.0));
}

#[test]
fn find_2_peaks_aspartic_acid() {
    let test = find_test("Aspartic acid");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(6),
                intensity_threshold: Some(300.0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);

    assert_eq!(peaks.len(), 2);

    // RT
    assert!(approx_eq(peaks[0].rt, 0.53, 0.05));
    assert!(approx_eq(peaks[1].rt, 0.88, 0.05));
    // Intensity
    assert!(approx_eq(peaks[0].intensity, 504.00, 100.0));
    assert!(approx_eq(peaks[1].intensity, 16462.00, 100.0));
    // Integral
    // assert!(approx_eq(peaks[0].integral, 12.1, 2.0));
    // assert!(approx_eq(peaks[1].integral, 546.6, 10.0));
}

#[test]
// #[ignore = "TODO"]
fn find_1_peak_proline_is() {
    let test = find_test("Proline IS");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(8),
                intensity_threshold: Some(400.0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 1);

    // RT
    assert!(approx_eq(peaks[0].rt, 2.479, 0.03));
    // // Intensity
    // assert!(approx_eq(peaks[0].intensity, 234.0, 100.0));
    // // Integral
    // assert!(approx_eq(peaks[0].integral, 4.6, 1.0));
}

#[test]
// #[ignore = "TODO"]
fn find_1_peak_glutamic_acid() {
    let test = find_test("Glutamic acid");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(8),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 1);

    // RT
    assert!(approx_eq(peaks[0].rt, 1.00, 0.03));
    // // Intensity
    // assert!(approx_eq(peaks[0].intensity, 234.0, 100.0));
    // // Integral
    // assert!(approx_eq(peaks[0].integral, 4.6, 1.0));
}

#[test]
fn find_1_peak_valine_is() {
    let test = find_test("Valine IS");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(8),
                intensity_threshold: Some(400.0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 1);

    // RT
    assert!(approx_eq(peaks[0].rt, 3.48, 0.03));
    // // Intensity
    // assert!(approx_eq(peaks[0].intensity, 234.0, 100.0));
    // // Integral
    // assert!(approx_eq(peaks[0].integral, 4.6, 1.0));
}

#[test]
fn find_4_peaks_sarcosine() {
    let test = find_test("Sarcosine");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(5),
                intensity_threshold: Some(400.0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 3);

    // RT
    // assert!(approx_eq(peaks[0].rt, 0.915, 0.03));
    assert!(approx_eq(peaks[0].rt, 2.083, 0.03));
    assert!(approx_eq(peaks[1].rt, 2.184, 0.03));
    assert!(approx_eq(peaks[2].rt, 2.385, 0.03));
    // Intensity
    // assert!(approx_eq(peaks[0].intensity, 384.0, 10.0));
    assert!(approx_eq(peaks[0].intensity, 708.0, 10.0));
    assert!(approx_eq(peaks[1].intensity, 1304.0, 10.0));
    assert!(approx_eq(peaks[2].intensity, 248550.0, 100.0));
}

#[test]
fn find_4_peaks_ce_18_1() {
    let test = find_test("CE(18:1)");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(5),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 2);

    // RT
    assert!(approx_eq(peaks[0].rt, 12.21, 0.03));
    assert!(approx_eq(peaks[1].rt, 12.47, 0.03));
    // Intensity
    assert!(approx_eq(peaks[0].intensity, 427518.0, 10.0));
    assert!(approx_eq(peaks[1].intensity, 296199.0, 10.0));
}

#[test]
fn find_4_peaks_ce_18_11() {
    let test = find_test("CE(18:11)");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(false),
                auto_baseline: Some(true),
                width_threshold: Some(15),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 1);

    // RT
    assert!(approx_eq(peaks[0].rt, 6.72, 0.03));
    // Intensity
    assert!(approx_eq(peaks[0].intensity, 20498768.0, 10.0));
}

#[test]
fn find_4_peaks_ce_18_12() {
    let test = find_test("CE(18:12)");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(15),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 1);

    // RT
    assert!(approx_eq(peaks[0].rt, 3.87, 0.03));
    // Intensity
    assert!(approx_eq(peaks[0].intensity, 73226.0, 10.0));
}

#[test]
#[ignore = "TODO: Includes band"]
fn find_6_peaks_285_114() {
    let test = find_test("COVp23_COV02260_285.114");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(8),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 7);

    assert!(approx_eq(peaks[0].rt, 3.487495, 0.03));
    assert!(approx_eq(peaks[0].intensity, 14196.0, 10.0));

    assert!(approx_eq(peaks[1].rt, 4.665315, 0.03));
    assert!(approx_eq(peaks[1].intensity, 5050.0, 10.0));

    assert!(approx_eq(peaks[2].rt, 4.904619, 0.03));
    assert!(approx_eq(peaks[2].intensity, 7752.0, 10.0));

    assert!(approx_eq(peaks[3].rt, 5.569014, 0.03));
    assert!(approx_eq(peaks[3].intensity, 6550.0, 10.0));

    assert!(approx_eq(peaks[4].rt, 5.635113, 0.03));
    assert!(approx_eq(peaks[4].intensity, 14166.0, 10.0));

    assert!(approx_eq(peaks[5].rt, 5.745212, 0.03));
    assert!(approx_eq(peaks[5].intensity, 15554.0, 10.0));

    assert!(approx_eq(peaks[6].rt, 5.954115, 0.03));
    assert!(approx_eq(peaks[6].intensity, 302250.0, 10.0));
}

#[test]
fn find_6_peaks_287_130() {
    let test = find_test("COVp23_COV02260_287.130");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(10),
                intensity_threshold: Some(400.0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 4);

    assert!(approx_eq(peaks[0].rt, 3.488905, 0.03));
    assert!(approx_eq(peaks[0].intensity, 11358.0, 10.0));

    assert!(approx_eq(peaks[1].rt, 5.367657, 0.03));
    assert!(approx_eq(peaks[1].intensity, 145096.0, 10.0));

    assert!(approx_eq(peaks[2].rt, 5.560452, 0.03));
    assert!(approx_eq(peaks[2].intensity, 72846.0, 10.0));

    assert!(approx_eq(peaks[3].rt, 5.626567, 0.03));
    assert!(approx_eq(peaks[3].intensity, 114256.0, 10.0));
}

#[test]
fn find_6_peaks_345_167() {
    let test = find_test("COVp23_COV02260_345.167");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 2);

    // #000
    assert!(approx_eq(peaks[0].rt, 2.593948, 0.03));
    assert!(approx_eq(peaks[0].intensity, 9610.0, 10.0));

    // #001
    assert!(approx_eq(peaks[1].rt, 5.947509, 0.03));
    assert!(approx_eq(peaks[1].intensity, 17716.0, 10.0));
}

#[test]
fn find_6_peaks_344_226() {
    let test = find_test("COVp23_COV02260_344.226");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(10),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 2);

    assert!(approx_eq(peaks[0].rt, 3.01, 0.03));
    assert!(approx_eq(peaks[1].rt, 5.60, 0.03));
    assert!(approx_eq(peaks[0].intensity, 1850.0, 10.0));
    assert!(approx_eq(peaks[1].intensity, 15356.0, 10.0));
}

#[test]
fn find_6_peaks_303_110() {
    let test = find_test("COVp23_COV02260_303.110");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 2);

    assert!(approx_eq(peaks[0].rt, 1.74, 0.03));
    assert!(approx_eq(peaks[1].rt, 2.03, 0.03));
    assert!(approx_eq(peaks[0].intensity, 12562.0, 10.0));
    assert!(approx_eq(peaks[1].intensity, 2464.0, 10.0));
}

#[test]
#[ignore = "TODO: Includes band"]
fn find_6_peaks_188_082() {
    let test = find_test("COVp23_COV02260_188.082");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 1);

    // #000
    assert!(approx_eq(peaks[0].rt, 5.61, 0.03));
    assert!(approx_eq(peaks[0].intensity, 15356.0, 10.0));
}

#[test]
fn find_6_peaks_379_239() {
    let test = find_test("COVp23_COV02260_379.239");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(10),
                intensity_threshold: Some(400.0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 1);

    // #000
    assert!(approx_eq(peaks[0].rt, 6.24, 0.03));
    assert!(approx_eq(peaks[0].intensity, 2900.0, 10.0));
}

#[test]
fn find_3_peaks_260_103() {
    let test = find_test("COVp23_COV02260_260.103");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(9),
                intensity_threshold: Some(400.0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 4);

    // RT
    assert!(approx_eq(peaks[0].rt, 2.087, 0.03));
    assert!(approx_eq(peaks[1].rt, 2.190, 0.03));
    assert!(approx_eq(peaks[2].rt, 2.386, 0.03));
    // Intensity
    assert!(approx_eq(peaks[0].intensity, 1316.0, 10.0));
    assert!(approx_eq(peaks[1].intensity, 916.0, 10.0));
    assert!(approx_eq(peaks[2].intensity, 107066.0, 100.0));
}

#[test]
fn find_1_peak_324_135() {
    let test = find_test("COVp26_COV02458_324.135");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_baseline: Some(false),
                auto_noise: Some(true),
                width_threshold: Some(5),
                intensity_threshold: Some(400.0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);

    assert_eq!(peaks.len(), 2);

    // RT
    assert!(approx_eq(peaks[0].rt, 1.88, 0.05));
    assert!(approx_eq(peaks[1].rt, 4.15, 0.05));
    // Intensity
    assert!(approx_eq(peaks[0].intensity, 6898.00, 100.0));
    assert!(approx_eq(peaks[1].intensity, 802.00, 10.0));
    // Integral
    // assert!(approx_eq(peaks[0].integral, 6219.72, 100.0));
}

#[test]
fn find_covp23_cov02260_peaks() {
    let test = find_test("COVp23_COV02260_peaks");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_baseline: Some(true),
                auto_noise: Some(true),
                width_threshold: Some(5),
                intensity_threshold: Some(400.0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);

    assert_eq!(peaks.len(), 30);

    assert!(approx_eq(peaks[0].rt, 0.850317, 0.05));
    assert!(approx_eq(peaks[1].rt, 1.037200, 0.05));
    assert!(approx_eq(peaks[2].rt, 1.290900, 0.05));
    assert!(approx_eq(peaks[3].rt, 1.738100, 0.05));
    assert!(approx_eq(peaks[4].rt, 1.789550, 0.05));
    assert!(approx_eq(peaks[5].rt, 1.905917, 0.05));
    assert!(approx_eq(peaks[6].rt, 1.940083, 0.05));
    assert!(approx_eq(peaks[7].rt, 2.035067, 0.05));
    assert!(approx_eq(peaks[8].rt, 2.138517, 0.05));
    assert!(approx_eq(peaks[9].rt, 2.190183, 0.05));
    assert!(approx_eq(peaks[10].rt, 2.267517, 0.05));
    assert!(approx_eq(peaks[11].rt, 2.385833, 0.05));
    assert!(approx_eq(peaks[12].rt, 2.638600, 0.05));
    assert!(approx_eq(peaks[13].rt, 2.890283, 0.05));
    assert!(approx_eq(peaks[14].rt, 3.170633, 0.05));
    assert!(approx_eq(peaks[15].rt, 3.439167, 0.05));
    assert!(approx_eq(peaks[16].rt, 3.473900, 0.05));
    assert!(approx_eq(peaks[17].rt, 3.548017, 0.05));
    assert!(approx_eq(peaks[18].rt, 4.017767, 0.05));
    assert!(approx_eq(peaks[19].rt, 4.069617, 0.05));
    assert!(approx_eq(peaks[20].rt, 4.155000, 0.05));
    assert!(approx_eq(peaks[21].rt, 4.224017, 0.05));
    assert!(approx_eq(peaks[22].rt, 4.352517, 0.05));
    assert!(approx_eq(peaks[23].rt, 4.899517, 0.05));
    assert!(approx_eq(peaks[24].rt, 5.366050, 0.05));
    assert!(approx_eq(peaks[25].rt, 5.557400, 0.05));
    assert!(approx_eq(peaks[26].rt, 5.626533, 0.05));
    assert!(approx_eq(peaks[27].rt, 5.743550, 0.05));
    assert!(approx_eq(peaks[28].rt, 5.954167, 0.05));
    assert!(approx_eq(peaks[29].rt, 6.070183, 0.05));

    assert!(approx_eq(peaks[0].intensity, 5392.0, 100.0));
    assert!(approx_eq(peaks[1].intensity, 15732.0, 100.0));
    assert!(approx_eq(peaks[2].intensity, 4324.0, 100.0));
    assert!(approx_eq(peaks[3].intensity, 12898.0, 100.0));
    assert!(approx_eq(peaks[4].intensity, 28074.0, 100.0));
    assert!(approx_eq(peaks[5].intensity, 129448.0, 100.0));
    assert!(approx_eq(peaks[6].intensity, 88962.0, 100.0));
    assert!(approx_eq(peaks[7].intensity, 3342.0, 100.0));
    assert!(approx_eq(peaks[8].intensity, 4240.0, 100.0));
    assert!(approx_eq(peaks[9].intensity, 3232.0, 100.0));
    assert!(approx_eq(peaks[10].intensity, 27536.0, 100.0));
    assert!(approx_eq(peaks[11].intensity, 321232.0, 100.0));
    assert!(approx_eq(peaks[12].intensity, 64716.0, 100.0));
    assert!(approx_eq(peaks[13].intensity, 26788.0, 100.0));
    assert!(approx_eq(peaks[14].intensity, 39464.0, 100.0));
    assert!(approx_eq(peaks[15].intensity, 75180.0, 100.0));
    assert!(approx_eq(peaks[16].intensity, 6602.0, 100.0));
    assert!(approx_eq(peaks[17].intensity, 24752.0, 100.0));
    assert!(approx_eq(peaks[18].intensity, 19978.0, 100.0));
    assert!(approx_eq(peaks[19].intensity, 53348.0, 100.0));
    assert!(approx_eq(peaks[20].intensity, 109812.0, 100.0));
    assert!(approx_eq(peaks[21].intensity, 53978.0, 100.0));
    assert!(approx_eq(peaks[22].intensity, 68226.0, 100.0));
    assert!(approx_eq(peaks[23].intensity, 4894.0, 100.0));
    assert!(approx_eq(peaks[24].intensity, 35884.0, 100.0));
    assert!(approx_eq(peaks[25].intensity, 21636.0, 100.0));
    assert!(approx_eq(peaks[26].intensity, 30796.0, 100.0));
    assert!(approx_eq(peaks[27].intensity, 4384.0, 100.0));
    assert!(approx_eq(peaks[28].intensity, 84156.0, 100.0));
    assert!(approx_eq(peaks[29].intensity, 3890.0, 100.0));
}

#[test]
fn find_rocit20_plate_3_ltr_19() {
    let test = find_test("ROCIT20-PLATE_3-LTR_19");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_baseline: Some(true),
                auto_noise: Some(false),
                width_threshold: Some(10),
                intensity_threshold: Some(400.0),
                ..Default::default()
            }),
            baseline_options: Some(BaselineOptions {
                baseline_window: Some(50.0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);

    assert_eq!(peaks.len(), 1);

    // RT
    assert!(approx_eq(peaks[0].rt, 6.92, 0.05));
    // Intensity
    assert!(approx_eq(peaks[0].intensity, 4317377.00, 100.0));
    // Integral
    // assert!(approx_eq(peaks[0].integral, 6219.72, 100.0));
}

#[test]
fn find_rocit20_plate_3_ltr_19_898_8_577_5() {
    let test = find_test("ROCIT20-PLATE_3-LTR_19_Q1=898.8 Q3=577.5");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_baseline: Some(true),
                auto_noise: Some(false),
                width_threshold: Some(5),
                intensity_threshold: Some(400.0),
                ..Default::default()
            }),
            baseline_options: Some(BaselineOptions {
                baseline_window: Some(50.0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);

    assert_eq!(peaks.len(), 1);

    // RT
    assert!(approx_eq(peaks[0].rt, 11.64, 0.05));
    // Intensity
    assert!(approx_eq(peaks[0].intensity, 35815576.00, 100.0));
    // Integral
    // assert!(approx_eq(peaks[0].integral, 6219.72, 100.0));
}

#[test]
fn find_rocit20_plate_3_rocit10263_692_6_369_4() {
    let test = find_test("PLATE_3-ROCIT10263_Q1=692.6 Q3=369.4");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_baseline: Some(true),
                auto_noise: Some(false),
                width_threshold: Some(8),
                intensity_threshold: Some(400.0),
                ..Default::default()
            }),
            baseline_options: Some(BaselineOptions {
                baseline_window: Some(50.0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);

    assert_eq!(peaks.len(), 1);

    // RT
    assert!(approx_eq(peaks[0].rt, 11.85, 0.05));
    // Intensity
    assert!(approx_eq(peaks[0].intensity, 8571918.00, 100.0));
    // Integral
    // assert!(approx_eq(peaks[0].integral, 6219.72, 100.0));
}

#[test]
fn find_rocit20_plate_3_rocit10059_952_8_679_5() {
    let test = find_test("PLATE_3-ROCIT10059_Q1=952.8 Q3=679.5");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_baseline: Some(true),
                auto_noise: Some(false),
                width_threshold: Some(4),
                intensity_threshold: Some(400.0),
                ..Default::default()
            }),
            baseline_options: Some(BaselineOptions {
                baseline_window: Some(50.0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);

    assert_eq!(peaks.len(), 2);

    // RT
    assert!(approx_eq(peaks[0].rt, 12.07, 0.05));
    // Intensity
    assert!(approx_eq(peaks[0].intensity, 376134.00, 100.0));
    // Integral
    // assert!(approx_eq(peaks[0].integral, 6219.72, 100.0));
}

#[test]
fn find_rocit20_plate_3_rocit10192_886_8_565_5() {
    let test = find_test("PLATE_3-ROCIT10192_Q1=886.8 Q3=565.5");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_baseline: Some(true),
                auto_noise: Some(false),
                width_threshold: Some(4),
                intensity_threshold: Some(400.0),
                ..Default::default()
            }),
            baseline_options: Some(BaselineOptions {
                baseline_window: Some(50.0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);

    assert_eq!(peaks.len(), 1);

    // RT
    assert!(approx_eq(peaks[0].rt, 11.41, 0.05));
    // Intensity
    assert!(approx_eq(peaks[0].intensity, 131972.00, 100.0));
    // Integral
    // assert!(approx_eq(peaks[0].integral, 6219.72, 100.0));
}

#[test]
fn find_rocit20_plate_3_ltr_19_902_8_565_5() {
    let test = find_test("PLATE_3-ROCIT10192_Q1=886.8 Q3=565.5");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_baseline: Some(true),
                auto_noise: Some(false),
                width_threshold: Some(4),
                intensity_threshold: Some(400.0),
                ..Default::default()
            }),
            baseline_options: Some(BaselineOptions {
                baseline_window: Some(50.0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);

    assert_eq!(peaks.len(), 1);

    // RT
    assert!(approx_eq(peaks[0].rt, 11.41, 0.05));
    // Intensity
    assert!(approx_eq(peaks[0].intensity, 131972.00, 100.0));
    // Integral
    // assert!(approx_eq(peaks[0].integral, 6219.72, 100.0));
}

#[test]
fn find_rocit20_plate_3_ltr_19_902_8_603_5() {
    let test = find_test("PLATE_3-LTR_19_Q1=902.8 Q3=603.5");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_baseline: Some(true),
                auto_noise: Some(false),
                width_threshold: Some(4),
                intensity_threshold: Some(400.0),
                ..Default::default()
            }),
            baseline_options: Some(BaselineOptions {
                baseline_window: Some(50.0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);

    assert_eq!(peaks.len(), 2);

    // RT
    assert!(approx_eq(peaks[0].rt, 12.22, 0.05));
    assert!(approx_eq(peaks[1].rt, 12.43, 0.05));
    // Intensity
    assert!(approx_eq(peaks[0].intensity, 59466468.00, 100.0));
    assert!(approx_eq(peaks[1].intensity, 60333712.00, 100.0));
    // Integral
    // assert!(approx_eq(peaks[0].integral, 6219.72, 100.0));
}

#[test]
#[ignore = "reason"]
fn find_rocit20_plate_3_10213_924_8_627_5() {
    let test = find_test("PLATE_3-rocit10213_Q1=924.8 Q3=627.5");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_baseline: Some(false),
                auto_noise: Some(false),
                width_threshold: Some(4),
                intensity_threshold: Some(400.0),
                ..Default::default()
            }),
            baseline_options: Some(BaselineOptions {
                baseline_window: Some(50.0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);

    assert_eq!(peaks.len(), 1);

    // RT
    assert!(approx_eq(peaks[0].rt, 11.59, 0.05));
    // Intensity
    assert!(approx_eq(peaks[0].intensity, 1329327.00, 100.0));
    // Integral
    // assert!(approx_eq(peaks[0].integral, 6219.72, 100.0));
}

#[test]
fn find_rocit20_plate_4_ltr_19_894_8_597_5() {
    let test = find_test("PLATE_4-ltr_19_Q1=894.8 Q3=597.5");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_baseline: Some(true),
                auto_noise: Some(false),
                width_threshold: Some(8),
                intensity_threshold: Some(400.0),
                ..Default::default()
            }),
            baseline_options: Some(BaselineOptions {
                baseline_window: Some(50.0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);

    assert_eq!(peaks.len(), 2);

    // RT
    assert!(approx_eq(peaks[0].rt, 10.39, 0.05));
    assert!(approx_eq(peaks[1].rt, 10.66, 0.05));
    // Intensity
    assert!(approx_eq(peaks[0].intensity, 173527.00, 100.0));
    assert!(approx_eq(peaks[1].intensity, 16561350.00, 100.0));
    // Integral
    // assert!(approx_eq(peaks[0].integral, 6219.72, 100.0));
}

#[test]
#[ignore = "reason"]
fn find_mauritius_cov2989_876_88_577_5() {
    let test = find_test("COVp298-COV19727_Q1=876.88 Q3=577.5");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_baseline: Some(true),
                auto_noise: Some(false),
                width_threshold: Some(4),
                intensity_threshold: Some(400.0),
                ..Default::default()
            }),
            baseline_options: Some(BaselineOptions {
                baseline_window: Some(50.0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);

    assert_eq!(peaks.len(), 3);

    // RT
    assert!(approx_eq(peaks[0].rt, 11.91, 0.05));
    assert!(approx_eq(peaks[1].rt, 12.19, 0.05));
    assert!(approx_eq(peaks[2].rt, 12.33, 0.05));
    // Intensity
    assert!(approx_eq(peaks[0].intensity, 1327232.00, 100.0));
    assert!(approx_eq(peaks[1].intensity, 111304408.00, 100.0));
    assert!(approx_eq(peaks[2].intensity, 110410112.00, 100.0));
    // Integral
    // assert!(approx_eq(peaks[0].integral, 6219.72, 100.0));
}

#[test]
fn find_mauritius_cov2989_815_7_184_1() {
    let test = find_test("COVp298-COV19727_Q1=815.7 Q3=184.1");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_baseline: Some(true),
                auto_noise: Some(false),
                width_threshold: Some(4),
                intensity_threshold: Some(400.0),
                ..Default::default()
            }),
            baseline_options: Some(BaselineOptions {
                baseline_window: Some(50.0),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);

    assert_eq!(peaks.len(), 3);

    // RT
    assert!(approx_eq(peaks[0].rt, 9.17, 0.05));
    assert!(approx_eq(peaks[1].rt, 9.35, 0.05));
    assert!(approx_eq(peaks[2].rt, 9.64, 0.05));
    // Intensity
    assert!(approx_eq(peaks[0].intensity, 18719300.00, 100.0));
    assert!(approx_eq(peaks[1].intensity, 31665968.00, 100.0));
    assert!(approx_eq(peaks[2].intensity, 116608344.00, 100.0));
}

#[test]
fn find_1_peaks_covp26_cov02499_324_1353() {
    let test = find_test("COVp26_COV02499_324.1353");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(10),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 3);

    assert!(approx_eq(peaks[0].rt, 1.90, 0.03));
    assert!(approx_eq(peaks[1].rt, 4.07, 0.03));
    assert!(approx_eq(peaks[2].rt, 4.16, 0.03));
    assert!(approx_eq(peaks[0].intensity, 12488.0, 10.0));
    assert!(approx_eq(peaks[1].intensity, 770.0, 10.0));
    assert!(approx_eq(peaks[2].intensity, 1372.0, 10.0));
}

#[test]
fn find_2_peaks_covp20_cov02037_302_1499() {
    let test = find_test("COVp20_COV02037_302.1499");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(10),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 2);

    assert!(approx_eq(peaks[0].rt, 4.07, 0.03));
    assert!(approx_eq(peaks[1].rt, 4.15, 0.03));
    assert!(approx_eq(peaks[0].intensity, 51404.0, 10.0));
    assert!(approx_eq(peaks[1].intensity, 107900.0, 10.0));
}

#[test]
fn find_1_peaks_covp24_cov02302_326_124() {
    let test = find_test("COVp24_COV02302_326.124");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(9),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 1);

    assert!(approx_eq(peaks[0].rt, 2.29, 0.03));
    assert!(approx_eq(peaks[0].intensity, 9226.0, 10.0));
}

#[test]
fn find_1_peaks_covp20_cov02022_274_1186() {
    let test = find_test("COVp20_COV02022_274.1186");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(3),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 2);

    assert!(approx_eq(peaks[0].rt, 2.54, 0.03));
    assert!(approx_eq(peaks[1].rt, 2.90, 0.03));
    assert!(approx_eq(peaks[0].intensity, 900.0, 10.0));
    assert!(approx_eq(peaks[1].intensity, 22336.0, 10.0));
}

#[test]
fn find_1_peaks_covp22_cal_8_6() {
    let test = find_test("COVp22_Cal_8_6");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(5),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 2);

    assert!(approx_eq(peaks[0].rt, 4.07, 0.03));
    assert!(approx_eq(peaks[1].rt, 4.15, 0.03));
    assert!(approx_eq(peaks[0].intensity, 3652.0, 10.0));
    assert!(approx_eq(peaks[1].intensity, 3824.0, 10.0));
}

#[test]
fn find_1_peaks_covp23_cal_4_120_237_1008() {
    let test = find_test("COVp23_Cal_4_120_237.1008");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(5),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 1);

    assert!(approx_eq(peaks[0].rt, 4.01, 0.03));
    assert!(approx_eq(peaks[0].intensity, 27586.0, 10.0));
}

#[test]
fn find_1_peaks_covp24_cal_2_12_274_1186() {
    let test = find_test("COVp24_Cal_2_12_274.1186");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(5),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    assert_eq!(peaks.len(), 3);

    assert!(approx_eq(peaks[0].rt, 2.46, 0.03));
    assert!(approx_eq(peaks[1].rt, 2.55, 0.03));
    assert!(approx_eq(peaks[2].rt, 2.88, 0.03));
    assert!(approx_eq(peaks[0].intensity, 183302.0, 10.0));
    assert!(approx_eq(peaks[1].intensity, 194902.0, 10.0));
    assert!(approx_eq(peaks[2].intensity, 195518.0, 10.0));
}

#[test]
fn find_1_peaks_covp26_cov02485_309_1671() {
    let test = find_test("COVp26_COV02485_309_1671");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(5),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 2);

    assert!(approx_eq(peaks[0].rt, 4.07, 0.03));
    assert!(approx_eq(peaks[1].rt, 4.16, 0.03));
    assert!(approx_eq(peaks[0].intensity, 15544.0, 10.0));
    assert!(approx_eq(peaks[1].intensity, 13270.0, 10.0));
}

#[test]
fn find_1_peaks_covp25_cal_7_117_288_1343() {
    let test = find_test("COVp25_Cal_7_117_288_1343");
    let peaks = find_peaks(
        &test.data,
        Some(FindPeaksOptions {
            filter_peaks_options: Some(FilterPeaksOptions {
                auto_noise: Some(true),
                auto_baseline: Some(true),
                width_threshold: Some(5),
                ..Default::default()
            }),
            ..Default::default()
        }),
    );

    dump_peaks(&peaks);
    assert_eq!(peaks.len(), 1);

    assert!(approx_eq(peaks[0].rt, 3.44, 0.03));
    assert!(approx_eq(peaks[0].intensity, 7932.0, 10.0));
}

// let peak: Peak = get_peak(
//     &test.data,
//     Roi {
//         rt: 2.55,
//         window: 0.05,
//     },
//     Some(FindPeaksOptions {
//         filter_peaks_options: Some(FilterPeaksOptions {
//             auto_noise: Some(true),
//             auto_baseline: Some(true),
//             width_threshold: Some(5),
//             ..Default::default()
//         }),
//         ..Default::default()
//     }),
// )
// .unwrap();
// let pk: Vec<Peak> = vec![peak];
// dump_peaks(&pk);
