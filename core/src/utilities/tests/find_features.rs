#[cfg(test)]
mod tests {
    use crate::utilities::{
        EicOptions,
        find_features::{
            Feature, FeaturePoint, FindFeaturesConfig, FindFeaturesOptions, MzTolerance,
            build_mz_grid, dedup_points, deduplicate_masses, group_by_mz, pick_best_per_rt_cluster,
            sort_features,
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
            mz_abs: v,
            ppm: 0.0,
        }
    }

    // --- build_mz_grid ---

    #[test]
    fn test_build_mz_grid_normal() {
        let grid = build_mz_grid(100.0, 102.0, 1.0);
        assert!(grid.len() >= 3);
        assert!(*grid.first().unwrap() >= 100.0);
        assert_eq!(*grid.last().unwrap(), 102.0);
    }

    #[test]
    fn test_build_mz_grid_nan_returns_empty() {
        assert!(build_mz_grid(f64::NAN, 100.0, 1.0).is_empty());
    }

    #[test]
    fn test_build_mz_grid_equal_bounds_returns_empty() {
        assert!(build_mz_grid(100.0, 100.0, 1.0).is_empty());
    }

    #[test]
    fn test_build_mz_grid_reversed_bounds_same_as_normal() {
        let forward = build_mz_grid(100.0, 102.0, 1.0);
        let reversed = build_mz_grid(102.0, 100.0, 1.0);
        assert_eq!(forward, reversed);
    }

    #[test]
    fn test_build_mz_grid_infinite_step_returns_endpoints() {
        let grid = build_mz_grid(100.0, 200.0, f64::INFINITY);
        assert_eq!(grid, vec![100.0, 200.0]);
    }

    #[test]
    fn test_build_mz_grid_last_point_is_exact_end() {
        let grid = build_mz_grid(100.0, 101.0, 0.5);
        assert_eq!(*grid.last().unwrap(), 101.0);
    }

    #[test]
    fn test_build_mz_grid_step_larger_than_range_returns_endpoints() {
        let grid = build_mz_grid(100.0, 101.0, 10.0);
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
            mz_abs: 0.001,
            ppm: 10.0,
        };
        let (lo, hi) = tol.window_at(1000.0);
        assert!((hi - lo) > 0.005);
    }

    #[test]
    fn test_mz_tolerance_abs_dominates_ppm_at_low_mz() {
        let tol = MzTolerance {
            mz_abs: 0.01,
            ppm: 1.0,
        };
        let (lo, hi) = tol.window_at(1.0);
        assert!((hi - lo - 0.02).abs() < 1e-10);
    }

    #[test]
    fn test_mz_tolerance_zero_ppm_uses_abs_only() {
        let tol = MzTolerance {
            mz_abs: 0.5,
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
        let sorted = deduplicate_masses(vec![100.0, 100.005, 200.0], opts.clone());
        let unsorted = deduplicate_masses(vec![100.005, 200.0, 100.0], opts);
        assert_eq!(sorted.len(), unsorted.len());
    }

    #[test]
    fn test_deduplicate_masses_two_elements_picks_upper_index() {
        let result = deduplicate_masses(vec![100.0, 100.005], eic_abs(0.01));
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], 100.005);
    }

    #[test]
    fn test_deduplicate_masses_three_elements_picks_middle() {
        let result = deduplicate_masses(vec![100.0, 100.005, 100.009], eic_abs(0.01));
        assert_eq!(result.len(), 1);
        assert_eq!(result[0], 100.005);
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

    // --- FindFeaturesConfig ---

    #[test]
    fn test_find_features_config_from_default_options() {
        let config = FindFeaturesConfig::from(FindFeaturesOptions::default());
        assert_eq!(config.scan_eic_options.ppm_tolerance, 10.0);
        assert_eq!(config.eic_options.ppm_tolerance, 20.0);
        assert_eq!(config.scan_width_threshold, 5);
        assert_eq!(config.mz_scan_grid.step_size, 0.006);
    }

    #[test]
    fn test_find_features_config_from_none_options_uses_defaults() {
        let config = FindFeaturesConfig::from(FindFeaturesOptions {
            scan_eic_options: None,
            eic_options: None,
            find_peaks: None,
            mz_scan_grid: None,
            scan_width_threshold: None,
        });
        assert_eq!(config.scan_eic_options.ppm_tolerance, 10.0);
        assert_eq!(config.scan_width_threshold, 5);
    }

    #[test]
    fn test_find_features_config_custom_values_passthrough() {
        let config = FindFeaturesConfig::from(FindFeaturesOptions {
            scan_eic_options: Some(EicOptions {
                ppm_tolerance: 5.0,
                mz_tolerance: 0.002,
                ..Default::default()
            }),
            eic_options: Some(EicOptions {
                ppm_tolerance: 15.0,
                mz_tolerance: 0.004,
                ..Default::default()
            }),
            scan_width_threshold: Some(3),
            find_peaks: None,
            mz_scan_grid: None,
        });
        assert_eq!(config.scan_eic_options.ppm_tolerance, 5.0);
        assert_eq!(config.eic_options.ppm_tolerance, 15.0);
        assert_eq!(config.scan_width_threshold, 3);
    }
}
