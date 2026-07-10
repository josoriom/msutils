#[cfg(test)]
mod tests {
    use crate::utilities::{
        EicOptions,
        calculate_eic::{lower_bound, mz_tolerance_for},
        find_features::{
            Feature, FeaturePoint, FindFeaturesOptions, MzTolerance, Scan, build_mz_grid,
            dedup_points, deduplicate_masses, group_by_mz, keep_masses_with_signal,
            max_eic_per_mass, pick_best_per_rt_cluster, sort_features,
        },
    };

    fn make_feature(mz: f64, rt: f64, intensity: f64) -> Feature {
        Feature {
            mz,
            rt,
            intensity,
            ..Default::default()
        }
    }

    fn make_fp(mz: f64, rt: f64, intensity: f64) -> FeaturePoint {
        FeaturePoint { rt, mz, intensity }
    }

    fn eic_abs(mz_tol: f64) -> EicOptions {
        EicOptions {
            ppm_tolerance: 0.0,
            mz_tolerance: mz_tol,
            ..Default::default()
        }
    }

    fn abs_tol(v: f64) -> MzTolerance {
        MzTolerance {
            mz_absolute: v,
            ppm: 0.0,
        }
    }

    // --- build_mz_grid ---

    #[test]
    fn test_build_mz_grid_normal() {
        let grid = build_mz_grid(100.0, 102.0, 1.0, 0.0);
        assert!(grid.len() >= 3);
        assert!(*grid.first().unwrap() >= 100.0);
        assert_eq!(*grid.last().unwrap(), 102.0);
    }

    #[test]
    fn test_build_mz_grid_nan_returns_empty() {
        assert!(build_mz_grid(f64::NAN, 100.0, 1.0, 0.0).is_empty());
    }

    #[test]
    fn test_build_mz_grid_equal_bounds_returns_empty() {
        assert!(build_mz_grid(100.0, 100.0, 1.0, 0.0).is_empty());
    }

    #[test]
    fn test_build_mz_grid_reversed_bounds_same_as_normal() {
        let forward = build_mz_grid(100.0, 102.0, 1.0, 0.0);
        let reversed = build_mz_grid(102.0, 100.0, 1.0, 0.0);
        assert_eq!(forward, reversed);
    }

    #[test]
    fn test_build_mz_grid_infinite_step_returns_endpoints() {
        let grid = build_mz_grid(100.0, 200.0, f64::INFINITY, 0.0);
        assert_eq!(grid, vec![100.0, 200.0]);
    }

    #[test]
    fn test_build_mz_grid_last_point_is_exact_end() {
        let grid = build_mz_grid(100.0, 101.0, 0.5, 0.0);
        assert_eq!(*grid.last().unwrap(), 101.0);
    }

    #[test]
    fn test_build_mz_grid_step_larger_than_range_returns_endpoints() {
        let grid = build_mz_grid(100.0, 101.0, 10.0, 0.0);
        assert_eq!(*grid.first().unwrap(), 100.0);
        assert_eq!(*grid.last().unwrap(), 101.0);
    }

    // --- MzTolerance ---

    #[test]
    fn test_mz_tolerance_are_close_within_abs() {
        assert!(abs_tol(0.01).are_close(100.0, 100.005));
    }

    #[test]
    fn test_mz_tolerance_are_close_outside_abs() {
        assert!(!abs_tol(0.01).are_close(100.0, 100.02));
    }

    #[test]
    fn test_mz_tolerance_are_close_to_ref_within() {
        assert!(abs_tol(0.01).are_close_to_ref(100.005, 100.0));
    }

    #[test]
    fn test_mz_tolerance_are_close_to_ref_outside() {
        assert!(!abs_tol(0.01).are_close_to_ref(100.02, 100.0));
    }

    #[test]
    fn test_mz_tolerance_are_close_to_ref_exact_match() {
        assert!(abs_tol(0.01).are_close_to_ref(100.0, 100.0));
    }

    #[test]
    fn test_mz_tolerance_window_at_symmetric() {
        let (lo, hi) = abs_tol(0.5).window_at(100.0);
        assert_eq!(lo, 99.5);
        assert_eq!(hi, 100.5);
    }

    #[test]
    fn test_mz_tolerance_ppm_dominates_abs_at_high_mz() {
        let tol = MzTolerance {
            mz_absolute: 0.001,
            ppm: 10.0,
        };
        let (lo, hi) = tol.window_at(1000.0);
        assert!((hi - lo) > 0.005);
    }

    #[test]
    fn test_mz_tolerance_abs_dominates_ppm_at_low_mz() {
        let tol = MzTolerance {
            mz_absolute: 0.01,
            ppm: 1.0,
        };
        let (lo, hi) = tol.window_at(1.0);
        assert!((hi - lo - 0.02).abs() < 1e-10);
    }

    #[test]
    fn test_mz_tolerance_zero_ppm_uses_abs_only() {
        let tol = MzTolerance {
            mz_absolute: 0.5,
            ppm: 0.0,
        };
        let (lo, hi) = tol.window_at(100.0);
        assert_eq!(lo, 99.5);
        assert_eq!(hi, 100.5);
    }

    // --- deduplicate_masses ---

    #[test]
    fn test_deduplicate_masses_empty() {
        assert!(deduplicate_masses(vec![], EicOptions::default()).is_empty());
    }

    #[test]
    fn test_deduplicate_masses_single() {
        assert_eq!(deduplicate_masses(vec![100.0], eic_abs(0.01)), vec![100.0]);
    }

    #[test]
    fn test_deduplicate_masses_no_duplicates() {
        let result = deduplicate_masses(vec![100.0, 200.0, 300.0], eic_abs(0.001));
        assert_eq!(result.len(), 3);
    }

    #[test]
    fn test_deduplicate_masses_merges_within_tolerance() {
        let result = deduplicate_masses(vec![100.0, 100.005, 200.0], eic_abs(0.01));
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_deduplicate_masses_unsorted_input_same_result() {
        let opts = eic_abs(0.01);
        let sorted = deduplicate_masses(vec![100.0, 100.005, 200.0], opts);
        let unsorted = deduplicate_masses(vec![100.005, 200.0, 100.0], opts);
        assert_eq!(sorted.len(), unsorted.len());
    }

    #[test]
    fn test_deduplicate_masses_two_elements_centers_on_mean() {
        let result = deduplicate_masses(vec![100.0, 100.005], eic_abs(0.01));
        assert_eq!(result.len(), 1);
        assert!((result[0] - 100.0025).abs() < 1e-9);
    }

    #[test]
    fn test_deduplicate_masses_cluster_centers_on_mean() {
        let result = deduplicate_masses(vec![100.0, 100.005, 100.009], eic_abs(0.01));
        assert_eq!(result.len(), 1);
        assert!((result[0] - (100.0 + 100.005 + 100.009) / 3.0).abs() < 1e-9);
    }

    // --- group_by_mz ---

    #[test]
    fn test_group_by_mz_empty() {
        let groups = group_by_mz(vec![], |a: f64, b: f64| (a - b).abs() <= 0.01);
        assert!(groups.is_empty());
    }

    #[test]
    fn test_group_by_mz_single() {
        let groups = group_by_mz(vec![make_fp(100.0, 1.0, 5.0)], |a, b| (a - b).abs() <= 0.01);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].len(), 1);
    }

    #[test]
    fn test_group_by_mz_all_same_mz() {
        let points = vec![make_fp(100.0, 1.0, 1.0), make_fp(100.0, 2.0, 2.0)];
        let groups = group_by_mz(points, |a, b| (a - b).abs() <= 0.01);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].len(), 2);
    }

    #[test]
    fn test_group_by_mz_two_groups() {
        let points = vec![
            make_fp(100.0, 1.0, 1.0),
            make_fp(100.005, 2.0, 2.0),
            make_fp(200.0, 3.0, 3.0),
        ];
        let groups = group_by_mz(points, |a, b| (a - b).abs() <= 0.01);
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].len(), 2);
        assert_eq!(groups[1].len(), 1);
    }

    #[test]
    fn test_group_by_mz_all_separate() {
        let points = vec![
            make_fp(100.0, 1.0, 1.0),
            make_fp(200.0, 2.0, 2.0),
            make_fp(300.0, 3.0, 3.0),
        ];
        let groups = group_by_mz(points, |a, b| (a - b).abs() <= 0.01);
        assert_eq!(groups.len(), 3);
    }

    // --- pick_best_per_rt_cluster ---

    #[test]
    fn test_pick_best_per_rt_cluster_empty() {
        assert!(pick_best_per_rt_cluster(vec![], 0.5).is_empty());
    }

    #[test]
    fn test_pick_best_per_rt_cluster_single() {
        let result = pick_best_per_rt_cluster(vec![make_fp(100.0, 1.0, 7.0)], 0.5);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].intensity, 7.0);
    }

    #[test]
    fn test_pick_best_per_rt_cluster_picks_max_intensity() {
        let group = vec![make_fp(100.0, 1.0, 5.0), make_fp(100.0, 1.1, 10.0)];
        let result = pick_best_per_rt_cluster(group, 0.5);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].intensity, 10.0);
    }

    #[test]
    fn test_pick_best_per_rt_cluster_two_separate_clusters() {
        let group = vec![make_fp(100.0, 1.0, 5.0), make_fp(100.0, 5.0, 3.0)];
        let result = pick_best_per_rt_cluster(group, 0.5);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_pick_best_per_rt_cluster_all_same_rt() {
        let group = vec![
            make_fp(100.0, 1.0, 5.0),
            make_fp(100.0, 1.0, 10.0),
            make_fp(100.0, 1.0, 3.0),
        ];
        let result = pick_best_per_rt_cluster(group, 0.5);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].intensity, 10.0);
    }

    #[test]
    fn test_pick_best_per_rt_cluster_all_separate() {
        let group = vec![
            make_fp(100.0, 1.0, 5.0),
            make_fp(100.0, 2.0, 3.0),
            make_fp(100.0, 3.0, 7.0),
        ];
        let result = pick_best_per_rt_cluster(group, 0.5);
        assert_eq!(result.len(), 3);
    }

    // --- dedup_points ---

    #[test]
    fn test_dedup_points_empty() {
        assert!(dedup_points(vec![], 0.01, 0.0, 0.5).is_empty());
    }

    #[test]
    fn test_dedup_points_single() {
        let result = dedup_points(vec![(1.0, 100.0, 5.0)], 0.01, 0.0, 0.5);
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn test_dedup_points_no_merging_needed() {
        let points = vec![(1.0, 100.0, 5.0), (2.0, 200.0, 3.0), (3.0, 300.0, 7.0)];
        assert_eq!(dedup_points(points, 0.001, 0.0, 0.1).len(), 3);
    }

    #[test]
    fn test_dedup_points_merges_nearby_rt_keeps_max_intensity() {
        let points = vec![(1.0, 100.0, 5.0), (1.05, 100.0, 10.0)];
        let result = dedup_points(points, 0.01, 0.0, 0.5);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].2, 10.0);
    }

    #[test]
    fn test_dedup_points_separate_mz_not_merged() {
        let points = vec![(1.0, 100.0, 5.0), (1.0, 200.0, 10.0)];
        assert_eq!(dedup_points(points, 0.01, 0.0, 0.5).len(), 2);
    }

    #[test]
    fn test_dedup_points_separate_rt_not_merged() {
        let points = vec![(1.0, 100.0, 5.0), (5.0, 100.0, 10.0)];
        assert_eq!(dedup_points(points, 0.01, 0.0, 0.5).len(), 2);
    }

    // --- sort_features ---

    #[test]
    fn test_sort_features_by_rt() {
        let features = vec![
            make_feature(100.0, 2.0, 500.0),
            make_feature(100.0, 1.0, 500.0),
        ];
        let result = sort_features(features);
        assert_eq!(result[0].rt, 1.0);
    }

    #[test]
    fn test_sort_features_same_rt_by_intensity_desc() {
        let features = vec![
            make_feature(100.0, 1.0, 100.0),
            make_feature(100.0, 1.0, 500.0),
        ];
        let result = sort_features(features);
        assert_eq!(result[0].intensity, 500.0);
    }

    #[test]
    fn test_sort_features_same_rt_same_intensity_by_mz_asc() {
        let features = vec![
            make_feature(200.0, 1.0, 500.0),
            make_feature(100.0, 1.0, 500.0),
        ];
        let result = sort_features(features);
        assert_eq!(result[0].mz, 100.0);
    }

    #[test]
    fn test_sort_features_empty() {
        assert!(sort_features(vec![]).is_empty());
    }

    #[test]
    fn test_sort_features_single() {
        let result = sort_features(vec![make_feature(100.0, 1.0, 500.0)]);
        assert_eq!(result.len(), 1);
    }

    // --- FindFeaturesOptions ---

    #[test]
    fn test_find_features_config_from_default_options() {
        let opts = FindFeaturesOptions::default();
        assert_eq!(opts.seed_eic_options.ppm_tolerance, 10.0);
        assert_eq!(opts.final_eic_options.ppm_tolerance, 20.0);
        assert_eq!(opts.min_seed_width_points, 5);
        assert_eq!(opts.mz_scan_grid.step, 0.006);
    }

    #[test]
    fn test_find_features_default_values() {
        let opts = FindFeaturesOptions::default();
        assert_eq!(opts.seed_eic_options.ppm_tolerance, 10.0);
        assert_eq!(opts.min_seed_width_points, 5);
    }

    #[test]
    fn test_find_features_custom_values_passthrough() {
        let opts = FindFeaturesOptions {
            seed_eic_options: EicOptions {
                ppm_tolerance: 5.0,
                mz_tolerance: 0.002,
                ..Default::default()
            },
            final_eic_options: EicOptions {
                ppm_tolerance: 15.0,
                mz_tolerance: 0.004,
                ..Default::default()
            },
            min_seed_width_points: 3,
            ..Default::default()
        };
        assert_eq!(opts.seed_eic_options.ppm_tolerance, 5.0);
        assert_eq!(opts.final_eic_options.ppm_tolerance, 15.0);
        assert_eq!(opts.min_seed_width_points, 3);
    }

    fn next_random(state: &mut u64) -> f64 {
        *state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*state >> 11) as f64) / ((1u64 << 53) as f64)
    }

    fn random_sorted(state: &mut u64, count: usize, low: f64, span: f64) -> Vec<f64> {
        let mut values: Vec<f64> = (0..count).map(|_| low + span * next_random(state)).collect();
        values.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        values
    }

    fn random_scan(state: &mut u64, points: usize) -> Scan {
        let mz = random_sorted(state, points, 100.0, 50.0);
        let intensity = (0..points).map(|_| next_random(state) * 1000.0).collect();
        Scan { mz, intensity }
    }

    fn reference_max_eic(scans: &[Scan], masses: &[f64], options: EicOptions) -> Vec<f64> {
        let lo: Vec<f64> = masses.iter().map(|&m| m - mz_tolerance_for(m, options)).collect();
        let hi: Vec<f64> = masses.iter().map(|&m| m + mz_tolerance_for(m, options)).collect();
        let mut max_value = vec![0.0f64; masses.len()];
        for scan in scans {
            let mz = &scan.mz;
            let intensity = &scan.intensity;
            let mut left = 0usize;
            let mut right = 0usize;
            for i in 0..masses.len() {
                while left < mz.len() && mz[left] < lo[i] {
                    left += 1;
                }
                if right < left {
                    right = left;
                }
                while right < mz.len() && mz[right] <= hi[i] {
                    right += 1;
                }
                let mut sum = 0.0;
                for k in left..right {
                    sum += intensity[k];
                }
                if sum > max_value[i] {
                    max_value[i] = sum;
                }
            }
        }
        max_value
    }

    fn reference_keep_masses(
        grid: &[f64],
        scans: &[Scan],
        min_intensity: f64,
        eic_options: EicOptions,
    ) -> Vec<f64> {
        let mut found_masses: Vec<f64> = Vec::new();
        for scan in scans {
            for (mass, intensity) in scan.mz.iter().zip(scan.intensity.iter()) {
                if *intensity >= min_intensity {
                    found_masses.push(*mass);
                }
            }
        }
        found_masses.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        grid.iter()
            .copied()
            .filter(|&mass| {
                let tolerance = mz_tolerance_for(mass, eic_options);
                let start = lower_bound(&found_masses, mass - tolerance);
                found_masses
                    .get(start)
                    .is_some_and(|&found| found <= mass + tolerance)
            })
            .collect()
    }

    fn random_count(state: &mut u64, max: f64) -> usize {
        1 + (next_random(state) * max) as usize
    }

    fn random_scans(state: &mut u64, count: usize) -> Vec<Scan> {
        (0..count)
            .map(|_| {
                let points = random_count(state, 40.0);
                random_scan(state, points)
            })
            .collect()
    }

    fn random_options(state: &mut u64) -> EicOptions {
        let ppm = next_random(state) * 40.0;
        let mz = next_random(state) * 0.02;
        EicOptions {
            ppm_tolerance: ppm,
            mz_tolerance: mz,
            ..Default::default()
        }
    }

    #[test]
    fn max_eic_per_mass_matches_reference_over_random_inputs() {
        let mut state = 0x1234_5678_9abc_def0u64;
        for _ in 0..300 {
            let scan_count = random_count(&mut state, 6.0);
            let scans = random_scans(&mut state, scan_count);
            let mass_count = random_count(&mut state, 30.0);
            let masses = random_sorted(&mut state, mass_count, 100.0, 50.0);
            let options = random_options(&mut state);
            let fast = max_eic_per_mass(&scans, &masses, options);
            let reference = reference_max_eic(&scans, &masses, options);
            assert_eq!(fast, reference);
        }
    }

    #[test]
    fn keep_masses_with_signal_matches_reference_over_random_inputs() {
        let mut state = 0x0fee_1dad_c0ff_ee11u64;
        for _ in 0..300 {
            let scan_count = random_count(&mut state, 6.0);
            let scans = random_scans(&mut state, scan_count);
            let grid_count = random_count(&mut state, 60.0);
            let grid = random_sorted(&mut state, grid_count, 100.0, 50.0);
            let min_intensity = next_random(&mut state) * 800.0;
            let options = random_options(&mut state);
            let fast = keep_masses_with_signal(&grid, &scans, min_intensity, options, 1);
            let reference = reference_keep_masses(&grid, &scans, min_intensity, options);
            assert_eq!(fast, reference);
        }
    }
}
