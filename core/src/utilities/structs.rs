use serde::{Deserialize, Serialize, Serializer};

pub fn ser_finite_f64<S: Serializer>(v: &f64, s: S) -> Result<S::Ok, S::Error> {
    s.serialize_f64(if v.is_finite() { *v } else { 0.0 })
}

#[derive(Debug, Deserialize)]
pub struct DataXY {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

#[derive(Clone, Copy, Debug)]
pub struct FromTo {
    pub from: f64,
    pub to: f64,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct Peak {
    #[serde(serialize_with = "ser_finite_f64")]
    pub from: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    pub to: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    pub rt: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    pub integral: f64,
    #[serde(serialize_with = "ser_finite_f64")]
    pub intensity: f64,
    pub n_points: usize,
    #[serde(serialize_with = "ser_finite_f64")]
    pub noise: f64,
    #[serde(skip)]
    pub r2: Option<f64>,
}

impl Default for Peak {
    fn default() -> Self {
        Self {
            from: 0.0,
            to: 0.0,
            rt: 0.0,
            integral: 0.0,
            intensity: 0.0,
            r2: None,
            n_points: 0,
            noise: 0.0,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Roi {
    pub rt: f64,
    pub range: f64,
}

impl Roi {
    pub fn new(rt: f64, range: f64) -> Self {
        Self { rt, range }
    }
}

#[derive(Clone, Debug)]
pub struct ChromRoi {
    pub id: String,
    pub sample_index: usize,
    pub rt: f64,
    pub range: f64,
}

impl ChromRoi {
    pub fn new(id: impl Into<String>, sample_index: usize, rt: f64, range: f64) -> Self {
        Self {
            id: id.into(),
            sample_index,
            rt,
            range,
        }
    }
}

#[derive(Clone, Debug, Deserialize)]
pub struct EicRoi {
    pub id: String,
    pub rt: f64,
    pub mz: f64,
    pub range: f64,
}

impl EicRoi {
    pub fn new(id: impl Into<String>, rt: f64, mz: f64, range: f64) -> Self {
        Self {
            id: id.into(),
            rt,
            mz,
            range,
        }
    }
}
