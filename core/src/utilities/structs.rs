use serde::Deserialize;

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

#[derive(Clone, Copy, Debug, Deserialize)]
pub struct Peak {
    pub from: f64,
    pub to: f64,
    pub rt: f64,
    pub integral: f64,
    pub intensity: f64,
    pub gaussian_r2: Option<f64>,
    pub n_points: usize,
    pub noise: f64,
}

impl Default for Peak {
    fn default() -> Self {
        Self {
            from: 0.0,
            to: 0.0,
            rt: 0.0,
            integral: 0.0,
            intensity: 0.0,
            gaussian_r2: None,
            n_points: 0,
            noise: 0.0,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Roi {
    pub rt: f64,
    pub half_width: f64,
}

impl Roi {
    pub fn new(rt: f64, half_width: f64) -> Self {
        Self { rt, half_width }
    }
}

#[derive(Clone, Debug)]
pub struct ChromRoi {
    pub id: String,
    pub sample_index: usize,
    pub rt: f64,
    pub half_width: f64,
}

impl ChromRoi {
    pub fn new(id: impl Into<String>, sample_index: usize, rt: f64, half_width: f64) -> Self {
        Self {
            id: id.into(),
            sample_index,
            rt,
            half_width,
        }
    }
}

#[derive(Clone, Debug, Deserialize)]
pub struct EicRoi {
    pub id: String,
    pub rt: f64,
    pub mz: f64,
    pub half_width: f64,
}

impl EicRoi {
    pub fn new(id: impl Into<String>, rt: f64, mz: f64, half_width: f64) -> Self {
        Self {
            id: id.into(),
            rt,
            mz,
            half_width,
        }
    }
}
