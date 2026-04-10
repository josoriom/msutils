struct FilterParams {
    num_targets: u32,
    num_scans:   u32,
    threshold:   f32,
}

@group(0) @binding(0) var<uniform>             params:           FilterParams;
@group(0) @binding(1) var<storage, read>       eic_matrix:       array<f32>;
@group(0) @binding(2) var<storage, read_write> survivor_count:   atomic<u32>;
@group(0) @binding(3) var<storage, read_write> survivor_indices: array<u32>;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
    let tidx = gid.x;
    if tidx >= params.num_targets { return; }

    let base = tidx * params.num_scans;
    var max_val = 0.0f;
    for (var s = 0u; s < params.num_scans; s++) {
        max_val = max(max_val, eic_matrix[base + s]);
    }
    if max_val <= params.threshold { return; }

    let slot = atomicAdd(&survivor_count, 1u);
    survivor_indices[slot] = tidx;
}