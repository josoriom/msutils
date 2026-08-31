use serde::Deserialize;


#[derive(Debug, Deserialize)]
pub struct DataXY {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

impl DataXY {
    pub fn empty() -> DataXY {
        Self {
            x: Vec::new(),
            y: Vec::new(),
        }
    }
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
    pub n_points: usize,
    pub noise: f64,
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
pub enum Roi {
    Peak {
        rt: f64,
        range: f64,
    },
    Eic {
        id: String,
        rt: f64,
        mz: f64,
        range: f64,
    },
    Chrom {
        id: String,
        sample_index: usize,
        rt: f64,
        range: f64,
    },
}

impl Roi {
    pub fn peak(rt: f64, range: f64) -> Self {
        Self::Peak { rt, range }
    }

    pub fn eic(id: impl Into<String>, rt: f64, mz: f64, range: f64) -> Self {
        Self::Eic {
            id: id.into(),
            rt,
            mz,
            range,
        }
    }

    pub fn chrom(id: impl Into<String>, sample_index: usize, rt: f64, range: f64) -> Self {
        Self::Chrom {
            id: id.into(),
            sample_index,
            rt,
            range,
        }
    }

    pub fn id(&self) -> &str {
        match self {
            Self::Peak { .. } => "",
            Self::Eic { id, .. } | Self::Chrom { id, .. } => id,
        }
    }
}
