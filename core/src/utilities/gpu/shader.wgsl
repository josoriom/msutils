struct Params {
    num_targets: u32,
    num_scans:   u32,
    ppm_tol:     f32,
    mz_tol:      f32,
}

@group(0) @binding(0) var<uniform>             params:  Params;
@group(0) @binding(1) var<storage, read>       scan_mz: array<f32>;
@group(0) @binding(2) var<storage, read>       scan_it: array<f32>;
@group(0) @binding(3) var<storage, read>       offsets: array<u32>;
@group(0) @binding(4) var<storage, read>       lengths: array<u32>;
@group(0) @binding(5) var<storage, read>       targets: array<f32>;
@group(0) @binding(6) var<storage, read_write> eic_out: array<f32>;

@compute @workgroup_size(64, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let tidx = gid.x;
    let sid  = gid.y;
    if tidx >= params.num_targets || sid >= params.num_scans {
        return;
    }

    let mz  = targets[tidx];
    let tol = max(params.ppm_tol * 1.0e-6 * mz, params.mz_tol);
    let lo  = mz - tol;
    let hi  = mz + tol;

    let base = offsets[sid];
    let len  = lengths[sid];
    var acc  = 0.0f;

    var left  = 0u;
    var right = len;
    while left < right {
        let mid = left + (right - left) / 2u;
        if scan_mz[base + mid] < lo {
            left = mid + 1u;
        } else {
            right = mid;
        }
    }

    var j = left;
    while j < len {
        let m = scan_mz[base + j];
        if m > hi { break; }
        acc += scan_it[base + j];
        j++;
    }

    eic_out[tidx * params.num_scans + sid] = acc;
}
