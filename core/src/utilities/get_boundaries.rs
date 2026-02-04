use crate::utilities::{closest_index, structs::DataXY};

#[derive(Clone, Copy, Debug)]
pub struct Boundary {
    pub index: Option<usize>,
    pub value: Option<f64>,
}

#[derive(Clone, Copy, Debug)]
pub struct Boundaries {
    pub from: Boundary,
    pub to: Boundary,
}

#[derive(Clone, Copy, Debug)]
pub struct BoundariesOptions {
    pub epsilon: f64,
    pub noise: f64,
    pub n_steps: usize,
    pub baseline_run: usize,
}

impl Default for BoundariesOptions {
    fn default() -> Self {
        Self {
            epsilon: 1e-5,
            noise: 0.0,
            n_steps: 3,
            baseline_run: 2,
        }
    }
}

pub fn get_boundaries(
    data: &DataXY,
    peak_x: f64,
    options: Option<BoundariesOptions>,
) -> Boundaries {
    let n = data.x.len();
    if n < 2 || n != data.y.len() {
        return Boundaries {
            from: Boundary {
                index: None,
                value: None,
            },
            to: Boundary {
                index: None,
                value: None,
            },
        };
    }

    let opts = options.unwrap_or_default();
    let idx = closest_index(&data.x, peak_x);

    let mut global_min = f64::INFINITY;
    for &vv in &data.y {
        let v = vv as f64;
        if v < global_min {
            global_min = v;
        }
    }

    let from = walk(
        &data.x,
        &data.y,
        idx,
        -1,
        opts.epsilon,
        opts.noise,
        opts.n_steps,
        opts.baseline_run,
        global_min,
    );

    let to = walk(
        &data.x,
        &data.y,
        idx,
        1,
        opts.epsilon,
        opts.noise,
        opts.n_steps,
        opts.baseline_run,
        global_min,
    );

    Boundaries { from, to }
}

fn walk(
    x: &[f64],
    y: &[f64],
    start: usize,
    direction: isize,
    epsilon: f64,
    noise_value: f64,
    n_steps: usize,
    baseline_run: usize,
    global_min: f64,
) -> Boundary {
    let n = x.len() as isize;
    if n < 2 {
        return Boundary {
            index: None,
            value: None,
        };
    }

    let dir = if direction > 0 { 1 } else { -1 };
    let noise = (global_min + epsilon).max(noise_value);

    let mut i = start as isize;

    let mut checking = false;
    let mut steps_up: usize = 0;
    let mut has_risen: bool = false;
    let mut valley_idx: isize = start as isize;
    let mut valley_val: f64 = y[start] as f64;
    let mut below_noise: bool = false;

    let mut below_noise_run_len: usize = 0;
    let mut below_noise_start: isize = 0;

    while i >= 0 && i + dir >= 0 && i + dir < n {
        let j = i + dir;
        let iu = i as usize;
        let ju = j as usize;

        let y_i = y[iu] as f64;
        let y_j = y[ju] as f64;

        if running_below_noise(
            y_j,
            noise,
            baseline_run,
            j,
            &mut below_noise_run_len,
            &mut below_noise_start,
        ) {
            let k = below_noise_start.clamp(0, n - 1) as usize;
            return Boundary {
                index: Some(k),
                value: Some(x[k]),
            };
        }

        let slope = compute_ratio_and_slope(x, y, iu, ju, dir, epsilon);

        if is_asc_or_flat(slope) {
            if !checking {
                checking = true;
                steps_up = 1;
                valley_idx = i;
                valley_val = y_i;
                below_noise = y_j <= noise;
            } else {
                steps_up += 1;
                if y_j <= noise {
                    below_noise = true;
                }
            }

            if steps_up >= n_steps {
                let end_val = y_j;
                let rise = end_val - valley_val;
                if !allow_rise(below_noise, end_val, noise, rise) {
                    let k = valley_idx.clamp(0, n - 1) as usize;
                    return Boundary {
                        index: Some(k),
                        value: Some(x[k]),
                    };
                }

                reset_state(
                    &mut checking,
                    &mut steps_up,
                    &mut has_risen,
                    &mut below_noise,
                );
            }
        } else {
            reset_state(
                &mut checking,
                &mut steps_up,
                &mut has_risen,
                &mut below_noise,
            );
        }

        i = j;
    }

    let edge = if dir > 0 { (n - 1) as usize } else { 0usize };
    Boundary {
        index: Some(edge),
        value: Some(x[edge]),
    }
}

fn running_below_noise(
    y_j: f64,
    noise: f64,
    baseline_run: usize,
    j: isize,
    below_noise_run_len: &mut usize,
    below_noise_start: &mut isize,
) -> bool {
    if y_j <= noise {
        if *below_noise_run_len == 0 {
            *below_noise_start = j;
        }
        *below_noise_run_len += 1;
        if *below_noise_run_len >= baseline_run {
            return true;
        }
    } else {
        *below_noise_run_len = 0;
    }
    false
}

fn reset_state(
    checking: &mut bool,
    steps_up: &mut usize,
    has_risen: &mut bool,
    below_noise: &mut bool,
) {
    *checking = false;
    *steps_up = 0;
    *has_risen = false;
    *below_noise = false;
}

fn is_asc_or_flat(slope: f64) -> bool {
    slope >= 0.0
}

#[inline]
fn allow_rise(below_noise: bool, end_val: f64, noise: f64, rise: f64) -> bool {
    !((below_noise && end_val <= noise) || (rise >= noise))
}

fn compute_ratio_and_slope(
    x: &[f64],
    y: &[f64],
    iu: usize,
    ju: usize,
    dir: isize,
    epsilon: f64,
) -> f64 {
    let dx = x[ju] - x[iu];
    let denom = if dx.abs() < epsilon {
        if dir > 0 { epsilon } else { -epsilon }
    } else {
        dx
    };
    let dy = (y[ju] as f64) - (y[iu] as f64);
    let slope_dir = (dy / denom) * if dir > 0 { 1.0 } else { -1.0 };
    slope_dir
}
