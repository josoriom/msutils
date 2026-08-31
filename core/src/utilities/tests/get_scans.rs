use std::sync::Arc;

use ionic::{
    WriteOptions,
    ion::{BytesSource, IonReader, ReadBytes, ReadOptions},
    mzml::structs::{
        BinaryDataArray, BinaryDataArrayList, CvParam, MzML, NumericArray, NumericType, Run, Scan,
        ScanList, Spectrum, SpectrumList,
    },
    write_mzml_to_ion,
};

use crate::{
    ffi::{FileSource, ParsedFile},
    utilities::{
        calculate_eic::{ScanQuery, TimeUnit, get_scans},
        structs::FromTo,
    },
};

const SCAN_MZ: [f64; 4] = [100.0, 499.999, 500.001, 900.0];
const SCAN_INTENSITY: [f64; 4] = [10.0, 1000.0, 900.0, 20.0];

fn cv_param(accession: &str, value: Option<&str>) -> CvParam {
    CvParam {
        accession: Some(accession.to_string()),
        name: accession.to_string(),
        value: value.map(str::to_string),
        ..Default::default()
    }
}

fn cv_param_with_unit(accession: &str, value: &str, unit_accession: &str) -> CvParam {
    CvParam {
        accession: Some(accession.to_string()),
        name: accession.to_string(),
        value: Some(value.to_string()),
        unit_accession: Some(unit_accession.to_string()),
        ..Default::default()
    }
}

fn binary_array(role_accession: &str, values: Vec<f64>) -> BinaryDataArray {
    BinaryDataArray {
        array_length: Some(values.len()),
        cv_params: vec![
            cv_param(role_accession, None),
            cv_param("MS:1000523", None),
            cv_param("MS:1000576", None),
        ],
        numeric_type: Some(NumericType::Float64),
        binary: Some(NumericArray::F64(values)),
        ..Default::default()
    }
}

fn spectrum(index: usize, rt: f64) -> Spectrum {
    Spectrum {
        id: format!("scan={}", index + 1),
        index: Some(index as u32),
        default_array_length: Some(SCAN_MZ.len()),
        ms_level: Some(1),
        cv_params: vec![cv_param("MS:1000511", Some("1"))],
        scan_list: Some(ScanList {
            count: Some(1),
            scans: vec![Scan {
                cv_params: vec![cv_param_with_unit(
                    "MS:1000016",
                    &rt.to_string(),
                    "UO:0000031",
                )],
                ..Default::default()
            }],
            ..Default::default()
        }),
        binary_data_array_list: Some(BinaryDataArrayList {
            count: Some(2),
            binary_data_arrays: vec![
                binary_array("MS:1000514", SCAN_MZ.to_vec()),
                binary_array("MS:1000515", SCAN_INTENSITY.to_vec()),
            ],
        }),
        ..Default::default()
    }
}

fn open_split_ion() -> ParsedFile {
    let mzml = MzML {
        run: Run {
            id: "test".to_string(),
            spectrum_list: Some(SpectrumList {
                count: Some(3),
                spectra: vec![spectrum(0, 1.0), spectrum(1, 2.0), spectrum(2, 3.0)],
                ..Default::default()
            }),
            ..Default::default()
        },
        ..Default::default()
    };

    let mut bytes = Vec::new();
    write_mzml_to_ion(
        &mzml,
        WriteOptions {
            compression_level: 0,
            force_f32: false,
            ..Default::default()
        },
        &mut bytes,
    )
    .expect("ion encode failed");
    let bytes = Arc::from(bytes.into_boxed_slice());
    let ion = IonReader::open_source(
        Arc::new(BytesSource::new(bytes)) as Arc<dyn ReadBytes>,
        ReadOptions::default(),
    )
    .expect("open ion failed");
    ParsedFile::new(FileSource::Lazy(Box::new(ion)))
}

#[test]
fn returns_every_segment_of_a_split_spectrum() {
    let mut file = open_split_ion();
    let (_, scans) = get_scans(
        &mut file,
        ScanQuery::RtRange(FromTo {
            from: 0.0,
            to: 10.0,
        }),
        TimeUnit::Minutes,
        1,
    );

    assert_eq!(scans.len(), 3);
    for scan in &scans {
        assert_eq!(scan.mz.as_ref(), SCAN_MZ.as_slice());
        assert_eq!(scan.intensity.as_ref(), SCAN_INTENSITY.as_slice());
    }
}

#[test]
fn rt_range_keeps_only_scans_in_window() {
    let mut file = open_split_ion();
    let (_, scans) = get_scans(
        &mut file,
        ScanQuery::RtRange(FromTo { from: 1.5, to: 2.5 }),
        TimeUnit::Minutes,
        1,
    );

    assert_eq!(scans.len(), 1);
    assert_eq!(scans[0].rt, 2.0);
    assert_eq!(scans[0].mz.as_ref(), SCAN_MZ.as_slice());
}

#[test]
fn closest_rt_returns_the_nearest_full_scan() {
    let mut file = open_split_ion();
    let (_, scans) = get_scans(&mut file, ScanQuery::ClosestRt(2.9), TimeUnit::Minutes, 1);

    assert_eq!(scans.len(), 1);
    assert_eq!(scans[0].rt, 3.0);
    assert_eq!(scans[0].mz.as_ref(), SCAN_MZ.as_slice());
}

#[test]
fn ms_level_filter_excludes_other_levels() {
    let mut file = open_split_ion();
    let (_, scans) = get_scans(
        &mut file,
        ScanQuery::RtRange(FromTo {
            from: 0.0,
            to: 10.0,
        }),
        TimeUnit::Minutes,
        2,
    );

    assert!(scans.is_empty());
}
