#[cfg(test)]
mod tests {
    use std::{cmp::Ordering, fs, sync::Arc};

    use ionic::{WriteOptions, mzml::structs::*, write_mzml_to_ion};

    use crate::utilities::{
        calculate_eic::{CentroidScan, EicOptions, SpectrumKind, SpectrumSummary},
        find_features::{Feature, FindFeaturesOptions, MzScanGrid, MzTolerance},
        get_features::{
            AlignmentOptions, ConsensusFeature, FeatureClusterer, SearchBounds,
            TaggedFeature, aggregate_into_consensus, assign_best_per_sample, collect_filled_slots,
            compute_search_bounds, dedup, get_features, require_minimum_frequency, resolve_cluster,
            weighted_centroid_mz,
        },
        math::median,
        mz_estimator::MzEstimatorKind,
        structs::FromTo,
    };

    fn make_scan(rt: f64, mz: Vec<f64>, intensity: Vec<f64>) -> CentroidScan {
        CentroidScan {
            rt,
            mz: Arc::from(mz.as_slice()),
            intensity: Arc::from(intensity.as_slice()),
            metadata: SpectrumSummary::unknown(),
        }
    }

    fn abs_eic(mz_tol: f64) -> EicOptions {
        EicOptions {
            ppm_tolerance: 0.0,
            mz_tolerance: mz_tol,
            ..Default::default()
        }
    }

    fn make_consensus(mz: f64, rt: f64, intensity: f64, n_samples: usize) -> ConsensusFeature {
        ConsensusFeature {
            mz,
            rt,
            from: rt - 0.05,
            to: rt + 0.05,
            intensity,
            integral: intensity * 0.1,
            frequency: 0.0,
            n_samples,
        }
    }

    fn abs_tol(v: f64) -> MzTolerance {
        MzTolerance {
            mz_absolute: v,
            ppm: 0.0,
        }
    }

    fn make_feature(mz: f64, rt: f64, intensity: f64) -> Feature {
        Feature {
            mz,
            rt,
            intensity,
            ..Default::default()
        }
    }

    fn make_feature_full(
        mz: f64,
        rt: f64,
        intensity: f64,
        n_points: usize,
        integral: f64,
    ) -> Feature {
        Feature {
            mz,
            rt,
            intensity,
            n_points,
            integral,
            ..Default::default()
        }
    }

    fn make_tagged(sample_index: usize, mz: f64, rt: f64, intensity: f64) -> TaggedFeature {
        TaggedFeature {
            sample_index,
            feature: make_feature(mz, rt, intensity),
        }
    }

    fn tagged_window(
        sample_index: usize,
        mz: f64,
        rt: f64,
        intensity: f64,
        from: f64,
        to: f64,
    ) -> TaggedFeature {
        TaggedFeature {
            sample_index,
            feature: Feature {
                mz,
                rt,
                intensity,
                from,
                to,
                ..Default::default()
            },
        }
    }

    fn default_clusterer() -> FeatureClusterer {
        FeatureClusterer {
            tolerance: MzTolerance {
                mz_absolute: 0.01,
                ppm: 0.0,
            },
            rt_tolerance: 0.1,
        }
    }

    fn default_bounds() -> SearchBounds {
        SearchBounds {
            target_mz: 100.0,
            center_rt: 1.0,
            rt_from: 0.9,
            rt_to: 1.1,
        }
    }

    #[test]
    fn test_median_two_elements() {
        assert_eq!(median(&mut [1.0, 3.0]), 2.0);
    }

    #[test]
    fn test_median_odd() {
        assert_eq!(median(&mut [3.0, 1.0, 2.0]), 2.0);
    }

    #[test]
    fn test_median_even() {
        assert_eq!(median(&mut [1.0, 2.0, 3.0, 4.0]), 2.5);
    }

    #[test]
    fn test_median_single() {
        assert_eq!(median(&mut [42.0]), 42.0);
    }

    #[test]
    fn test_median_unsorted_input() {
        assert_eq!(median(&mut [5.0, 1.0, 3.0]), 3.0);
    }





    #[test]
    fn test_clusterer_empty_input() {
        let clusters = default_clusterer().cluster(vec![]);
        assert!(clusters.is_empty());
    }

    #[test]
    fn test_clusterer_single_item() {
        let clusters = default_clusterer().cluster(vec![make_tagged(0, 100.0, 1.0, 500.0)]);
        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0].len(), 1);
    }

    #[test]
    fn test_clusterer_groups_close_features() {
        let clusterer = default_clusterer();
        let items = vec![
            tagged_window(0, 100.0, 1.0, 500.0, 0.95, 1.10),
            tagged_window(1, 100.005, 1.05, 400.0, 1.00, 1.15),
            tagged_window(0, 200.0, 2.0, 300.0, 1.95, 2.05),
        ];
        let clusters = clusterer.cluster(items);
        assert_eq!(clusters.len(), 2);
        assert_eq!(clusters[0].len(), 2);
        assert_eq!(clusters[1].len(), 1);
    }

    #[test]
    fn test_clusterer_separates_distant_rt() {
        let clusterer = default_clusterer();
        let items = vec![
            make_tagged(0, 100.0, 1.0, 500.0),
            make_tagged(1, 100.005, 5.0, 400.0),
        ];
        let clusters = clusterer.cluster(items);
        assert_eq!(clusters.len(), 2);
    }

    #[test]
    fn test_clusterer_just_outside_mz_boundary_splits() {
        let items = vec![
            make_tagged(0, 100.0, 1.0, 500.0),
            make_tagged(1, 100.011, 1.0, 400.0),
        ];
        let clusters = default_clusterer().cluster(items);
        assert_eq!(clusters.len(), 2);
    }

    #[test]
    fn test_clusterer_mz_within_tolerance_groups() {
        let items = vec![
            tagged_window(0, 100.0, 1.0, 500.0, 0.95, 1.05),
            tagged_window(1, 100.005, 1.0, 400.0, 0.95, 1.05),
        ];
        let clusters = default_clusterer().cluster(items);
        assert_eq!(clusters.len(), 1);
    }

    #[test]
    fn test_clusterer_mz_outside_tolerance_splits() {
        let items = vec![
            make_tagged(0, 100.0, 1.0, 500.0),
            make_tagged(1, 100.02, 1.0, 400.0),
        ];
        let clusters = default_clusterer().cluster(items);
        assert_eq!(clusters.len(), 2);
    }

    #[test]
    fn test_clusterer_rt_within_tolerance_groups() {
        let items = vec![
            tagged_window(0, 100.0, 1.0, 500.0, 0.95, 1.10),
            tagged_window(1, 100.005, 1.05, 400.0, 1.00, 1.15),
        ];
        let clusters = default_clusterer().cluster(items);
        assert_eq!(clusters.len(), 1);
    }

    #[test]
    fn test_clusterer_rt_outside_tolerance_splits() {
        let items = vec![
            make_tagged(0, 100.0, 1.0, 500.0),
            make_tagged(1, 100.005, 1.2, 400.0),
        ];
        let clusters = default_clusterer().cluster(items);
        assert_eq!(clusters.len(), 2);
    }

    #[test]
    fn test_clusterer_just_outside_rt_boundary_splits() {
        let items = vec![
            make_tagged(0, 100.0, 1.0, 500.0),
            make_tagged(1, 100.005, 1.101, 400.0),
        ];
        let clusters = default_clusterer().cluster(items);
        assert_eq!(clusters.len(), 2);
    }

    #[test]
    fn test_assign_best_per_sample_picks_highest_intensity() {
        let cluster = vec![
            make_tagged(0, 100.0, 1.0, 100.0),
            make_tagged(0, 100.0, 1.0, 500.0),
            make_tagged(1, 100.0, 1.0, 200.0),
        ];
        let slots = assign_best_per_sample(cluster, 2);
        assert_eq!(slots[0].as_ref().unwrap().intensity, 500.0);
        assert_eq!(slots[1].as_ref().unwrap().intensity, 200.0);
    }

    #[test]
    fn test_assign_best_per_sample_empty_cluster() {
        let slots = assign_best_per_sample(vec![], 3);
        assert!(slots.iter().all(|s| s.is_none()));
    }

    #[test]
    fn test_assign_best_per_sample_single_sample() {
        let cluster = vec![make_tagged(0, 100.0, 1.0, 500.0)];
        let slots = assign_best_per_sample(cluster, 1);
        assert_eq!(slots[0].as_ref().unwrap().intensity, 500.0);
    }

    #[test]
    fn test_assign_best_per_sample_missing_sample_is_none() {
        let cluster = vec![make_tagged(0, 100.0, 1.0, 500.0)];
        let slots = assign_best_per_sample(cluster, 3);
        assert!(slots[1].is_none());
        assert!(slots[2].is_none());
    }

    #[test]
    fn test_compute_search_bounds_all_none_returns_none() {
        let slots: Vec<Option<Feature>> = vec![None, None];
        assert!(compute_search_bounds(&slots, 0.05).is_none());
    }

    #[test]
    fn test_compute_search_bounds_single_slot() {
        let slots = vec![Some(make_feature(100.0, 1.0, 500.0))];
        let bounds = compute_search_bounds(&slots, 0.05).unwrap();
        assert_eq!(bounds.target_mz, 100.0);
        assert_eq!(bounds.center_rt, 1.0);
        assert_eq!(bounds.rt_from, 0.95);
        assert_eq!(bounds.rt_to, 1.05);
    }

    #[test]
    fn test_compute_search_bounds_two_slots_upper_index() {
        let slots = vec![
            Some(make_feature(100.0, 1.0, 500.0)),
            Some(make_feature(102.0, 2.0, 300.0)),
        ];
        let bounds = compute_search_bounds(&slots, 0.05).unwrap();
        assert_eq!(bounds.target_mz, 102.0);
        assert_eq!(bounds.center_rt, 2.0);
    }

    #[test]
    fn test_compute_search_bounds_correct_rt_window() {
        let slots = vec![
            Some(make_feature(100.0, 1.0, 500.0)),
            Some(make_feature(100.0, 1.0, 400.0)),
            Some(make_feature(100.0, 1.0, 300.0)),
        ];
        let bounds = compute_search_bounds(&slots, 0.1).unwrap();
        assert!((bounds.rt_from - 0.9).abs() < 1e-10);
        assert!((bounds.rt_to - 1.1).abs() < 1e-10);
    }

    #[test]
    fn test_require_minimum_frequency_exact_threshold() {
        let hits = vec![make_feature(100.0, 1.0, 500.0)];
        assert!(require_minimum_frequency(hits, 1).is_some());
    }

    #[test]
    fn test_require_minimum_frequency_passes() {
        let hits = vec![
            make_feature(100.0, 1.0, 500.0),
            make_feature(100.0, 1.0, 400.0),
        ];
        assert!(require_minimum_frequency(hits, 2).is_some());
    }

    #[test]
    fn test_require_minimum_frequency_fails() {
        let hits = vec![make_feature(100.0, 1.0, 500.0)];
        assert!(require_minimum_frequency(hits, 3).is_none());
    }

    #[test]
    fn test_require_minimum_frequency_empty() {
        assert!(require_minimum_frequency(vec![], 1).is_none());
    }

    #[test]
    fn test_collect_filled_slots_filters_nones() {
        let slots = vec![
            Some(make_feature(100.0, 1.0, 500.0)),
            None,
            Some(make_feature(200.0, 2.0, 300.0)),
        ];
        assert_eq!(collect_filled_slots(slots).len(), 2);
    }

    #[test]
    fn test_collect_filled_slots_all_none() {
        let slots: Vec<Option<Feature>> = vec![None, None, None];
        assert!(collect_filled_slots(slots).is_empty());
    }

    #[test]
    fn test_collect_filled_slots_all_some() {
        let slots = vec![
            Some(make_feature(100.0, 1.0, 500.0)),
            Some(make_feature(200.0, 2.0, 300.0)),
        ];
        assert_eq!(collect_filled_slots(slots).len(), 2);
    }

    #[test]
    fn test_aggregate_frequency() {
        let hits = vec![
            make_feature(100.0, 1.0, 500.0),
            make_feature(100.0, 1.0, 400.0),
            make_feature(100.0, 1.0, 300.0),
        ];
        let result = aggregate_into_consensus(hits, &default_bounds(), 3);
        assert_eq!(result.n_samples, 3);
    }

    #[test]
    fn test_aggregate_from_to_from_bounds() {
        let hits = vec![make_feature(100.0, 1.0, 500.0)];
        let bounds = SearchBounds {
            target_mz: 100.0,
            center_rt: 1.0,
            rt_from: 0.8,
            rt_to: 1.2,
        };
        let result = aggregate_into_consensus(hits, &bounds, 1);
        assert_eq!(result.from, 0.8);
        assert_eq!(result.to, 1.2);
    }

    #[test]
    fn test_aggregate_intensity_even_median() {
        let hits = vec![
            make_feature(100.0, 1.0, 300.0),
            make_feature(100.0, 1.0, 500.0),
        ];
        let result = aggregate_into_consensus(hits, &default_bounds(), 3);
        assert_eq!(result.intensity, 400.0);
    }

    #[test]
    fn test_aggregate_mz_odd_median() {
        let hits = vec![
            make_feature(99.0, 1.0, 500.0),
            make_feature(100.0, 1.0, 400.0),
            make_feature(101.0, 1.0, 300.0),
        ];
        let result = aggregate_into_consensus(hits, &default_bounds(), 3);
        assert_eq!(result.mz, 100.0);
    }

    #[test]
    fn test_aggregate_frequency_ratio() {
        // 3 hits out of 6 total samples → frequency = 0.5
        let hits = vec![
            make_feature_full(100.0, 1.0, 500.0, 2, 10.0),
            make_feature_full(100.0, 1.0, 400.0, 4, 20.0),
            make_feature_full(100.0, 1.0, 300.0, 6, 30.0),
        ];
        let result = aggregate_into_consensus(hits, &default_bounds(), 6);
        assert!((result.frequency - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_aggregate_integral_median() {
        let hits = vec![
            make_feature_full(100.0, 1.0, 500.0, 2, 10.0),
            make_feature_full(100.0, 1.0, 400.0, 4, 20.0),
            make_feature_full(100.0, 1.0, 300.0, 6, 30.0),
        ];
        let result = aggregate_into_consensus(hits, &default_bounds(), 3);
        assert!((result.integral - 20.0).abs() < 1e-10);
    }

    #[test]
    fn test_aggregate_correct_frequency() {
        let hits = vec![
            make_feature(100.0, 1.0, 500.0),
            make_feature(100.01, 1.05, 400.0),
            make_feature(99.99, 0.95, 600.0),
        ];
        let result = aggregate_into_consensus(hits, &default_bounds(), 3);
        assert_eq!(result.n_samples, 3);
        assert!((result.mz - 100.0).abs() < 0.01);
    }

    // --- weighted_centroid_mz ---

    #[test]
    fn test_wcmz_peak_at_target_returns_target() {
        let scans = vec![make_scan(5.0, vec![100.0], vec![1000.0])];
        let times = vec![5.0];
        let result = weighted_centroid_mz(&scans, &times, 100.0, abs_eic(0.01), 4.9, 5.1);
        assert!((result.unwrap() - 100.0).abs() < 1e-9);
    }

    #[test]
    fn test_wcmz_peak_offset_low_returns_actual_centroid() {
        // peak at 99.995, target at 100.0, window covers it → should return 99.995, not 100.0
        let scans = vec![make_scan(5.0, vec![99.995], vec![1000.0])];
        let times = vec![5.0];
        let result = weighted_centroid_mz(&scans, &times, 100.0, abs_eic(0.01), 4.9, 5.1);
        assert!((result.unwrap() - 99.995).abs() < 1e-9);
    }

    #[test]
    fn test_wcmz_peak_offset_high_returns_actual_centroid() {
        // peak at 100.005, target at 100.0 → should return 100.005
        let scans = vec![make_scan(5.0, vec![100.005], vec![1000.0])];
        let times = vec![5.0];
        let result = weighted_centroid_mz(&scans, &times, 100.0, abs_eic(0.01), 4.9, 5.1);
        assert!((result.unwrap() - 100.005).abs() < 1e-9);
    }

    #[test]
    fn test_wcmz_heavy_left_pulls_centroid_below_midpoint() {
        // 90% of intensity on low side → centroid should be left of 100.0
        let scans = vec![make_scan(5.0, vec![99.995, 100.005], vec![9000.0, 1000.0])];
        let times = vec![5.0];
        let result = weighted_centroid_mz(&scans, &times, 100.0, abs_eic(0.01), 4.9, 5.1).unwrap();
        assert!(result < 100.0, "centroid {result} should be below 100.0");
        assert!((result - 99.996).abs() < 1e-6);
    }

    #[test]
    fn test_wcmz_heavy_right_pulls_centroid_above_midpoint() {
        // 90% of intensity on high side → centroid should be right of 100.0
        let scans = vec![make_scan(5.0, vec![99.995, 100.005], vec![1000.0, 9000.0])];
        let times = vec![5.0];
        let result = weighted_centroid_mz(&scans, &times, 100.0, abs_eic(0.01), 4.9, 5.1).unwrap();
        assert!(result > 100.0, "centroid {result} should be above 100.0");
        assert!((result - 100.004).abs() < 1e-6);
    }

    #[test]
    fn test_wcmz_peak_outside_mz_window_returns_none() {
        // peak at 100.02, window is ±0.01 → miss
        let scans = vec![make_scan(5.0, vec![100.02], vec![1000.0])];
        let times = vec![5.0];
        assert!(weighted_centroid_mz(&scans, &times, 100.0, abs_eic(0.01), 4.9, 5.1).is_none());
    }

    #[test]
    fn test_wcmz_scan_outside_rt_window_returns_none() {
        // scan at rt=10.0, window [4.9, 5.1] → no signal
        let scans = vec![make_scan(10.0, vec![100.0], vec![1000.0])];
        let times = vec![10.0];
        assert!(weighted_centroid_mz(&scans, &times, 100.0, abs_eic(0.01), 4.9, 5.1).is_none());
    }

    #[test]
    fn test_wcmz_zero_target_returns_none() {
        let scans = vec![make_scan(5.0, vec![0.0], vec![1000.0])];
        let times = vec![5.0];
        assert!(weighted_centroid_mz(&scans, &times, 0.0, abs_eic(0.01), 4.9, 5.1).is_none());
    }

    #[test]
    fn test_wcmz_nan_target_returns_none() {
        let scans = vec![make_scan(5.0, vec![100.0], vec![1000.0])];
        let times = vec![5.0];
        assert!(weighted_centroid_mz(&scans, &times, f64::NAN, abs_eic(0.01), 4.9, 5.1).is_none());
    }

    #[test]
    fn test_wcmz_multi_scan_weighted_across_rt() {
        // two scans at different rt, both within window — centroid is intensity-weighted
        let scans = vec![
            make_scan(4.95, vec![99.996], vec![8000.0]),
            make_scan(5.05, vec![100.004], vec![2000.0]),
        ];
        let times = vec![4.95, 5.05];
        let result = weighted_centroid_mz(&scans, &times, 100.0, abs_eic(0.01), 4.9, 5.1).unwrap();
        // expected: (99.996*8000 + 100.004*2000) / 10000 = 99.9976
        let expected = (99.996 * 8000.0 + 100.004 * 2000.0) / 10000.0;
        assert!((result - expected).abs() < 1e-9);
    }

    // --- dedup ---

    #[test]
    fn test_dedup_empty_input() {
        let result = dedup(vec![], &abs_tol(0.005), 0.05);
        assert!(result.is_empty());
    }

    #[test]
    fn test_dedup_single_feature_kept() {
        let input = vec![make_consensus(100.0, 5.0, 1000.0, 5)];
        let result = dedup(input, &abs_tol(0.005), 0.05);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].n_samples, 5);
    }

    #[test]
    fn test_dedup_higher_frequency_wins_large() {
        // 9 vs 6 — the bug scenario
        let high = make_consensus(100.000, 5.0, 7000.0, 9);
        let low = make_consensus(100.002, 5.0, 7100.0, 6);
        let result = dedup(vec![high, low], &abs_tol(0.005), 0.05);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].n_samples, 9);
    }

    #[test]
    fn test_dedup_higher_frequency_wins_small() {
        // 2 vs 1 — boundary case
        let high = make_consensus(200.000, 3.0, 500.0, 2);
        let low = make_consensus(200.001, 3.0, 600.0, 1);
        let result = dedup(vec![high, low], &abs_tol(0.005), 0.05);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].n_samples, 2);
    }

    #[test]
    fn test_dedup_intensity_tiebreaker_same_frequency() {
        let strong = make_consensus(100.000, 5.0, 9000.0, 4);
        let weak = make_consensus(100.002, 5.0, 3000.0, 4);
        let result = dedup(vec![weak, strong], &abs_tol(0.005), 0.05);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].intensity as u64, 9000);
    }

    #[test]
    fn test_dedup_intensity_tiebreaker_reversed_order() {
        // same as above but inserted in different order
        let strong = make_consensus(100.000, 5.0, 9000.0, 4);
        let weak = make_consensus(100.002, 5.0, 3000.0, 4);
        let result = dedup(vec![strong, weak], &abs_tol(0.005), 0.05);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].intensity as u64, 9000);
    }

    #[test]
    fn test_dedup_chain_three_within_tolerance_best_survives() {
        // three features all within tolerance of each other
        let a = make_consensus(100.000, 5.0, 5000.0, 8);
        let b = make_consensus(100.002, 5.0, 4000.0, 5);
        let c = make_consensus(100.004, 5.0, 3000.0, 3);
        let result = dedup(vec![a, b, c], &abs_tol(0.005), 0.05);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].n_samples, 8);
    }

    #[test]
    fn test_dedup_just_outside_mz_tolerance_both_kept() {
        // 0.006 > abs_tol(0.005) → separate features
        let a = make_consensus(100.000, 5.0, 5000.0, 5);
        let b = make_consensus(100.006, 5.0, 4000.0, 4);
        let result = dedup(vec![a, b], &abs_tol(0.005), 0.05);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_dedup_just_outside_rt_tolerance_both_kept() {
        let a = make_consensus(100.000, 5.000, 5000.0, 5);
        let b = make_consensus(100.001, 5.051, 4000.0, 4); // rt diff = 0.051 > 0.05
        let result = dedup(vec![a, b], &abs_tol(0.005), 0.05);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_dedup_just_inside_mz_tolerance_lower_freq_removed() {
        // 0.004 < abs_tol(0.005) → same feature
        let a = make_consensus(100.000, 5.0, 5000.0, 7);
        let b = make_consensus(100.004, 5.0, 6000.0, 3);
        let result = dedup(vec![a, b], &abs_tol(0.005), 0.05);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].n_samples, 7);
    }

    #[test]
    fn test_dedup_just_inside_rt_tolerance_lower_freq_removed() {
        let a = make_consensus(100.000, 5.000, 5000.0, 6);
        let b = make_consensus(100.001, 5.049, 6000.0, 2); // rt diff = 0.049 < 0.05
        let result = dedup(vec![a, b], &abs_tol(0.005), 0.05);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].n_samples, 6);
    }

    fn gaussian(rt: f64, mu: f64, sigma: f64, amp: f64) -> f64 {
        amp * (-0.5 * ((rt - mu) / sigma).powi(2)).exp()
    }

    fn cv(accession: &str, value: Option<&str>) -> CvParam {
        CvParam {
            accession: Some(accession.to_string()),
            name: accession.to_string(),
            value: value.map(str::to_string),
            ..Default::default()
        }
    }

    fn cv_unit(accession: &str, value: &str, unit_accession: &str) -> CvParam {
        CvParam {
            accession: Some(accession.to_string()),
            name: accession.to_string(),
            value: Some(value.to_string()),
            unit_accession: Some(unit_accession.to_string()),
            ..Default::default()
        }
    }

    fn binary_array(role_accession: &str, values: Vec<f64>) -> BinaryDataArray {
        let len = values.len();
        BinaryDataArray {
            array_length: Some(len),
            cv_params: vec![
                cv(role_accession, None),
                cv("MS:1000523", None), // 64-bit float
                cv("MS:1000576", None), // no compression
            ],
            numeric_type: Some(NumericType::Float64),
            binary: Some(NumericArray::F64(values)),
            ..Default::default()
        }
    }

    fn make_spectrum(index: usize, rt_min: f64, mz: Vec<f64>, intensity: Vec<f64>) -> Spectrum {
        let len = mz.len();
        Spectrum {
            id: format!("scan={}", index + 1),
            index: Some(index as u32),
            default_array_length: Some(len),
            ms_level: Some(1),
            cv_params: vec![cv("MS:1000511", Some("1"))], // ms level
            scan_list: Some(ScanList {
                count: Some(1),
                scans: vec![Scan {
                    cv_params: vec![cv_unit(
                        "MS:1000016", // scan start time
                        &rt_min.to_string(),
                        "UO:0000031", // minute
                    )],
                    ..Default::default()
                }],
                ..Default::default()
            }),
            binary_data_array_list: Some(BinaryDataArrayList {
                count: Some(2),
                binary_data_arrays: vec![
                    binary_array("MS:1000514", mz),
                    binary_array("MS:1000515", intensity),
                ],
            }),
            ..Default::default()
        }
    }

    fn build_mzml_single_peak(
        peak_mz: f64,
        peak_rt: f64,
        peak_amp: f64,
        rt_min: f64,
        rt_max: f64,
        n_scans: usize,
    ) -> MzML {
        let sigma = 0.05;
        let mut spectra = Vec::with_capacity(n_scans);
        for i in 0..n_scans {
            let rt = rt_min + (rt_max - rt_min) * i as f64 / (n_scans - 1) as f64;
            let intensity = gaussian(rt, peak_rt, sigma, peak_amp).max(0.0);
            spectra.push(make_spectrum(i, rt, vec![peak_mz], vec![intensity]));
        }
        MzML {
            run: Run {
                id: "test".to_string(),
                spectrum_list: Some(SpectrumList {
                    count: Some(spectra.len()),
                    spectra,
                    ..Default::default()
                }),
                ..Default::default()
            },
            ..Default::default()
        }
    }

    fn build_mzml_two_peaks(
        mz_a: f64,
        mz_b: f64,
        peak_rt: f64,
        amp: f64,
        rt_min: f64,
        rt_max: f64,
        n_scans: usize,
    ) -> MzML {
        let sigma = 0.05;
        let mut spectra = Vec::with_capacity(n_scans);
        for i in 0..n_scans {
            let rt = rt_min + (rt_max - rt_min) * i as f64 / (n_scans - 1) as f64;
            let ia = gaussian(rt, peak_rt, sigma, amp).max(0.0);
            let ib = gaussian(rt, peak_rt, sigma, amp * 0.8).max(0.0);
            spectra.push(make_spectrum(i, rt, vec![mz_a, mz_b], vec![ia, ib]));
        }
        MzML {
            run: Run {
                id: "test".to_string(),
                spectrum_list: Some(SpectrumList {
                    count: Some(spectra.len()),
                    spectra,
                    ..Default::default()
                }),
                ..Default::default()
            },
            ..Default::default()
        }
    }

    fn build_mzml_masses(
        mzs: &[f64],
        peak_rt: f64,
        amp: f64,
        rt_min: f64,
        rt_max: f64,
        n_scans: usize,
    ) -> MzML {
        let sigma = 0.05;
        let mut spectra = Vec::with_capacity(n_scans);
        for i in 0..n_scans {
            let rt = rt_min + (rt_max - rt_min) * i as f64 / (n_scans - 1) as f64;
            let height = gaussian(rt, peak_rt, sigma, amp).max(0.0);
            let intensities: Vec<f64> = (0..mzs.len())
                .map(|k| height * (1.0 - 0.15 * k as f64))
                .collect();
            spectra.push(make_spectrum(i, rt, mzs.to_vec(), intensities));
        }
        MzML {
            run: Run {
                id: "test".to_string(),
                spectrum_list: Some(SpectrumList {
                    count: Some(spectra.len()),
                    spectra,
                    ..Default::default()
                }),
                ..Default::default()
            },
            ..Default::default()
        }
    }

    fn build_mzml_profile_peak(
        center_mz: f64,
        mz_step: f64,
        n_mz: usize,
        mz_sigma: f64,
        peak_rt: f64,
        amp: f64,
        rt_min: f64,
        rt_max: f64,
        n_scans: usize,
    ) -> MzML {
        let rt_sigma = 0.05;
        let half = (n_mz / 2) as f64;
        let mut spectra = Vec::with_capacity(n_scans);
        for i in 0..n_scans {
            let rt = rt_min + (rt_max - rt_min) * i as f64 / (n_scans - 1) as f64;
            let rt_height = gaussian(rt, peak_rt, rt_sigma, amp).max(0.0);
            let mut mzs = Vec::with_capacity(n_mz);
            let mut intensities = Vec::with_capacity(n_mz);
            for k in 0..n_mz {
                let m = center_mz + (k as f64 - half) * mz_step;
                mzs.push(m);
                intensities.push(rt_height * gaussian(m, center_mz, mz_sigma, 1.0));
            }
            let mut spectrum = make_spectrum(i, rt, mzs, intensities);
            spectrum.cv_params.push(cv("MS:1000128", None));
            spectra.push(spectrum);
        }
        MzML {
            run: Run {
                id: "test".to_string(),
                spectrum_list: Some(SpectrumList {
                    count: Some(spectra.len()),
                    spectra,
                    ..Default::default()
                }),
                ..Default::default()
            },
            ..Default::default()
        }
    }

    fn build_mzml_peak_and_flat(
        peak_mz: f64,
        flat_mz: f64,
        peak_rt: f64,
        amp: f64,
        rt_min: f64,
        rt_max: f64,
        n_scans: usize,
    ) -> MzML {
        let sigma = 0.05;
        let flat_level = amp * 0.3;
        let mut spectra = Vec::with_capacity(n_scans);
        for i in 0..n_scans {
            let rt = rt_min + (rt_max - rt_min) * i as f64 / (n_scans - 1) as f64;
            let height = gaussian(rt, peak_rt, sigma, amp).max(0.0);
            spectra.push(make_spectrum(
                i,
                rt,
                vec![peak_mz, flat_mz],
                vec![height, flat_level],
            ));
        }
        MzML {
            run: Run {
                id: "test".to_string(),
                spectrum_list: Some(SpectrumList {
                    count: Some(spectra.len()),
                    spectra,
                    ..Default::default()
                }),
                ..Default::default()
            },
            ..Default::default()
        }
    }

    fn to_ion_bytes(mzml: &MzML) -> Vec<u8> {
        let mut bytes = Vec::new();
        write_mzml_to_ion(
            mzml,
            WriteOptions {
                compression_level: 0,
                force_f32: false,
                ..Default::default()
            },
            &mut bytes,
        )
        .expect("ion encode should succeed");
        bytes
    }

    fn write_ion_dir(tag: &str, samples: &[Vec<u8>]) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("quantion_test_{tag}"));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).expect("create temp dir");
        for (i, bytes) in samples.iter().enumerate() {
            let path = dir.join(format!("sample_{i:02}.ion"));
            fs::write(&path, bytes).expect("write ion file");
        }
        dir
    }

    fn alignment_cfg() -> AlignmentOptions {
        AlignmentOptions {
            mz_tolerance: MzTolerance {
                mz_absolute: 0.005,
                ppm: 20.0,
            },
            rt_tolerance: 0.15,
            min_samples: 1,
            eic_options: EicOptions {
                ppm_tolerance: 20.0,
                mz_tolerance: 0.005,
                ..Default::default()
            },
            peak_options: None,
            mz_estimator: MzEstimatorKind::MedianMzApex,
            ..Default::default()
        }
    }

    fn feature_at(mz: f64, rt: f64, from: f64, to: f64) -> Feature {
        Feature {
            mz,
            rt,
            from,
            to,
            intensity: 0.0,
            integral: 0.0,
            n_points: 0,
            noise: 0.0,
        }
    }

    fn bounds_at(target_mz: f64, center_rt: f64, rt_from: f64, rt_to: f64) -> SearchBounds {
        SearchBounds {
            target_mz,
            center_rt,
            rt_from,
            rt_to,
        }
    }

    #[test]
    fn test_resolve_cluster_present_keeps_feature_and_reads_apex() {
        let scans = vec![
            (4.9, vec![100.001], vec![10.0]),
            (5.0, vec![100.002], vec![90.0]),
            (5.1, vec![100.003], vec![15.0]),
        ];
        let existing = feature_at(100.0, 5.0, 4.9, 5.1);
        let bounds = bounds_at(100.0, 5.0, 4.85, 5.15);
        let resolved = resolve_cluster(
            &scans,
            Some(&existing),
            &bounds,
            99.99,
            100.01,
            20.0,
            None,
            SpectrumKind::Centroid,
        );
        assert!(resolved.feature.is_none());
        let apex = resolved.apex.expect("apex measured");
        assert!((apex.mz - 100.002).abs() < 1e-9);
        assert!((apex.intensity - 90.0).abs() < 1e-9);
    }

    #[test]
    fn test_resolve_cluster_apex_uses_only_peak_window() {
        let scans = vec![
            (4.0, vec![100.050], vec![1000.0]),
            (5.0, vec![100.002], vec![90.0]),
            (6.0, vec![100.090], vec![1000.0]),
        ];
        let existing = feature_at(100.0, 5.0, 4.9, 5.1);
        let bounds = bounds_at(100.0, 5.0, 4.85, 5.15);
        let resolved = resolve_cluster(
            &scans,
            Some(&existing),
            &bounds,
            99.0,
            101.0,
            20.0,
            None,
            SpectrumKind::Centroid,
        );
        let apex = resolved.apex.expect("apex measured");
        assert!((apex.mz - 100.002).abs() < 1e-9);
    }

    #[test]
    fn test_resolve_cluster_missing_flat_signal_returns_nothing() {
        let scans = vec![
            (4.9, vec![100.002], vec![0.0]),
            (5.0, vec![100.002], vec![0.0]),
            (5.1, vec![100.002], vec![0.0]),
        ];
        let bounds = bounds_at(100.0, 5.0, 4.85, 5.15);
        let resolved = resolve_cluster(
            &scans,
            None,
            &bounds,
            99.99,
            100.01,
            20.0,
            None,
            SpectrumKind::Centroid,
        );
        assert!(resolved.feature.is_none());
        assert!(resolved.apex.is_none());
    }

    #[test]
    fn test_resolve_cluster_missing_with_peak_fills_and_measures_apex() {
        let peak_mz = 300.0;
        let scans: Vec<(f64, Vec<f64>, Vec<f64>)> = (0..41)
            .map(|i| {
                let rt = 4.0 + 2.0 * i as f64 / 40.0;
                (rt, vec![peak_mz], vec![gaussian(rt, 5.0, 0.1, 1000.0)])
            })
            .collect();
        let bounds = bounds_at(peak_mz, 5.0, 4.7, 5.3);
        let resolved = resolve_cluster(
            &scans,
            None,
            &bounds,
            299.99,
            300.01,
            20.0,
            None,
            SpectrumKind::Centroid,
        );

        let feature = resolved
            .feature
            .expect("a missing sample with a real peak should be filled");
        assert!(
            feature.from < 5.0 && feature.to > 5.0,
            "filled peak [{:.3}, {:.3}] should bracket the apex rt 5.0",
            feature.from,
            feature.to
        );
        assert!(
            (feature.mz - peak_mz).abs() < 1e-9,
            "filled feature mz should keep the slot target"
        );
        let apex = resolved.apex.expect("apex measured for the filled peak");
        assert!((apex.mz - peak_mz).abs() < 1e-9);
    }

    fn feature_cfg_for_mz(mz: f64) -> FindFeaturesOptions {
        FindFeaturesOptions {
            mz_scan_grid: MzScanGrid {
                min_mz: (mz - 2.0).max(40.0),
                max_mz: (mz + 2.0).min(1000.0),
                ..Default::default()
            },
            ..Default::default()
        }
    }

    fn feature_cfg_for_mz_range(low: f64, high: f64) -> FindFeaturesOptions {
        FindFeaturesOptions {
            mz_scan_grid: MzScanGrid {
                min_mz: (low - 2.0).max(40.0),
                max_mz: (high + 2.0).min(1000.0),
                ..Default::default()
            },
            ..Default::default()
        }
    }

    #[test]
    fn test_full_cluster_no_fills_detects_single_feature() {
        let n_samples = 5;
        let peak_mz = 200.0;
        let peak_rt = 5.0;

        let samples: Vec<Vec<u8>> = (0..n_samples)
            .map(|_| {
                let mzml = build_mzml_single_peak(peak_mz, peak_rt, 10_000.0, 4.0, 6.0, 60);
                to_ion_bytes(&mzml)
            })
            .collect();

        let dir = write_ion_dir("full_cluster", &samples);
        let features = get_features(
            dir.to_str().unwrap(),
            FromTo { from: 3.0, to: 7.0 },
            feature_cfg_for_mz(peak_mz),
            alignment_cfg(),
            1,
        )
        .expect("get_features should succeed");

        let _ = fs::remove_dir_all(&dir);

        assert!(
            !features.is_empty(),
            "should detect at least one feature, got none"
        );
        let best = features.iter().max_by_key(|f| f.n_samples).unwrap();
        assert!(
            (best.mz - peak_mz).abs() < 0.01,
            "consensus mz {:.4} should be near {peak_mz}",
            best.mz
        );
        assert!(
            (best.rt - peak_rt).abs() < 0.2,
            "consensus rt {:.3} should be near {peak_rt}",
            best.rt
        );
        assert_eq!(
            best.n_samples, n_samples,
            "all {n_samples} samples should be in the cluster"
        );
    }

    #[test]
    fn test_full_cluster_high_mz_detects_feature() {
        let n_samples = 4;
        let peak_mz = 993.674;
        let peak_rt = 9.0;

        let samples: Vec<Vec<u8>> = (0..n_samples)
            .map(|_| {
                let mzml = build_mzml_single_peak(peak_mz, peak_rt, 8_000.0, 8.0, 10.0, 60);
                to_ion_bytes(&mzml)
            })
            .collect();

        let dir = write_ion_dir("full_cluster_high_mz", &samples);
        let features = get_features(
            dir.to_str().unwrap(),
            FromTo {
                from: 7.0,
                to: 11.0,
            },
            feature_cfg_for_mz(peak_mz),
            alignment_cfg(),
            1,
        )
        .expect("get_features should succeed");

        let _ = fs::remove_dir_all(&dir);

        assert!(!features.is_empty(), "should detect at least one feature");
        let best = features.iter().max_by_key(|f| f.n_samples).unwrap();
        assert!(
            (best.mz - peak_mz).abs() < 0.02,
            "consensus mz {:.4} should be near {peak_mz:.4}",
            best.mz
        );
    }

    #[test]
    fn test_filled_mz_matches_measured_centroid() {
        let n_real = 8;
        let peak_mz = 300.0;
        let peak_rt = 5.0;

        let samples: Vec<Vec<u8>> = (0..n_real)
            .map(|_| {
                let mzml = build_mzml_single_peak(peak_mz, peak_rt, 12_000.0, 4.0, 6.0, 60);
                to_ion_bytes(&mzml)
            })
            .collect();

        let dir = write_ion_dir("fill_centroid", &samples);
        let features = get_features(
            dir.to_str().unwrap(),
            FromTo { from: 3.5, to: 6.5 },
            feature_cfg_for_mz(peak_mz),
            alignment_cfg(),
            1,
        )
        .expect("get_features should succeed");

        let _ = fs::remove_dir_all(&dir);

        assert!(!features.is_empty());
        let best = features.iter().max_by_key(|f| f.n_samples).unwrap();
        assert!(
            (best.mz - peak_mz).abs() < 0.01,
            "filled consensus mz {:.4} should be near the real peak mz {peak_mz:.4}",
            best.mz
        );
    }

    #[test]
    fn test_filled_mz_matches_measured_centroid_high_mz() {
        let peak_mz = 750.123;
        let peak_rt = 7.5;

        let samples: Vec<Vec<u8>> = (0..6)
            .map(|_| {
                let mzml = build_mzml_single_peak(peak_mz, peak_rt, 9_000.0, 6.5, 8.5, 60);
                to_ion_bytes(&mzml)
            })
            .collect();

        let dir = write_ion_dir("fill_centroid_high_mz", &samples);
        let features = get_features(
            dir.to_str().unwrap(),
            FromTo { from: 6.0, to: 9.0 },
            feature_cfg_for_mz(peak_mz),
            alignment_cfg(),
            1,
        )
        .expect("get_features should succeed");

        let _ = fs::remove_dir_all(&dir);

        assert!(!features.is_empty());
        let best = features.iter().max_by_key(|f| f.n_samples).unwrap();
        assert!(
            (best.mz - peak_mz).abs() < 0.02,
            "consensus mz {:.4} should be near {peak_mz:.4}",
            best.mz
        );
    }

    #[test]
    fn test_frequency_threshold_filters_rare_feature() {
        let common_mz = 150.0;
        let rare_mz = 155.0;
        let peak_rt = 5.0;
        let n_samples = 5;

        let samples: Vec<Vec<u8>> = (0..n_samples)
            .enumerate()
            .map(|(i, _)| {
                let mzml = if i == 0 {
                    // first sample has both peaks
                    build_mzml_two_peaks(common_mz, rare_mz, peak_rt, 10_000.0, 4.0, 6.0, 60)
                } else {
                    // the other samples have only the common peak
                    build_mzml_single_peak(common_mz, peak_rt, 10_000.0, 4.0, 6.0, 60)
                };
                to_ion_bytes(&mzml)
            })
            .collect();

        let dir = write_ion_dir("freq_threshold", &samples);
        let mut cfg = alignment_cfg();
        cfg.min_samples = 3; // require presence in at least 3 samples

        let features = get_features(
            dir.to_str().unwrap(),
            FromTo { from: 3.0, to: 7.0 },
            feature_cfg_for_mz_range(common_mz, rare_mz),
            cfg,
            1,
        )
        .expect("get_features should succeed");

        let _ = fs::remove_dir_all(&dir);

        assert!(
            features.iter().any(|f| (f.mz - common_mz).abs() < 0.01),
            "common feature at mz≈{common_mz} should survive frequency filter"
        );
        assert!(
            !features
                .iter()
                .any(|f| (f.mz - rare_mz).abs() < 0.01 && f.n_samples < 3),
            "rare feature at mz≈{rare_mz} below threshold should be filtered out"
        );
    }

    #[test]
    fn test_frequency_threshold_exact_match_kept() {
        let peak_mz = 180.0;
        let peak_rt = 4.0;
        let n_samples = 4;

        let samples: Vec<Vec<u8>> = (0..n_samples)
            .map(|_| {
                let mzml = build_mzml_single_peak(peak_mz, peak_rt, 8_000.0, 3.0, 5.0, 60);
                to_ion_bytes(&mzml)
            })
            .collect();

        let dir = write_ion_dir("freq_exact", &samples);
        let mut cfg = alignment_cfg();
        cfg.min_samples = n_samples; // threshold = exact count

        let features = get_features(
            dir.to_str().unwrap(),
            FromTo { from: 2.5, to: 5.5 },
            feature_cfg_for_mz(peak_mz),
            cfg,
            1,
        )
        .expect("get_features should succeed");

        let _ = fs::remove_dir_all(&dir);

        assert!(
            features.iter().any(|f| (f.mz - peak_mz).abs() < 0.01),
            "feature present in all {n_samples} samples should pass threshold={n_samples}"
        );
    }

    #[test]
    fn test_two_disjoint_peaks_both_reported() {
        let mz_a = 200.000;
        let mz_b = 210.000; // 10 Da apart — no overlap
        let peak_rt = 5.0;
        let n_samples = 4;

        let samples: Vec<Vec<u8>> = (0..n_samples)
            .map(|_| {
                let mzml = build_mzml_two_peaks(mz_a, mz_b, peak_rt, 9_000.0, 4.0, 6.0, 60);
                to_ion_bytes(&mzml)
            })
            .collect();

        let dir = write_ion_dir("two_disjoint", &samples);
        let features = get_features(
            dir.to_str().unwrap(),
            FromTo { from: 3.0, to: 7.0 },
            feature_cfg_for_mz_range(mz_a, mz_b),
            alignment_cfg(),
            1,
        )
        .expect("get_features should succeed");

        let _ = fs::remove_dir_all(&dir);

        let has_a = features.iter().any(|f| (f.mz - mz_a).abs() < 0.05);
        let has_b = features.iter().any(|f| (f.mz - mz_b).abs() < 0.05);
        assert!(has_a, "compound A at mz≈{mz_a} should be reported");
        assert!(has_b, "compound B at mz≈{mz_b} should be reported");
    }

    #[test]
    fn test_two_close_peaks_outside_tolerance_both_kept() {
        let mz_a = 500.000;
        let mz_b = 500.030;
        let peak_rt = 6.0;
        let n_samples = 3;

        let samples: Vec<Vec<u8>> = (0..n_samples)
            .map(|_| {
                let mzml = build_mzml_two_peaks(mz_a, mz_b, peak_rt, 8_000.0, 5.0, 7.0, 60);
                to_ion_bytes(&mzml)
            })
            .collect();

        let dir = write_ion_dir("two_close", &samples);
        let features = get_features(
            dir.to_str().unwrap(),
            FromTo { from: 4.0, to: 8.0 },
            feature_cfg_for_mz_range(mz_a, mz_b),
            alignment_cfg(),
            1,
        )
        .expect("get_features should succeed");

        let _ = fs::remove_dir_all(&dir);

        let has_a = features.iter().any(|f| (f.mz - mz_a).abs() < 0.015);
        let has_b = features.iter().any(|f| (f.mz - mz_b).abs() < 0.015);
        assert!(has_a, "compound A at mz≈{mz_a} should be reported");
        assert!(has_b, "compound B at mz≈{mz_b} should be reported");
    }

    #[test]
    fn test_consensus_feature_keeps_exact_values() {
        let f = ConsensusFeature {
            mz: f64::NAN,
            rt: f64::INFINITY,
            from: f64::NEG_INFINITY,
            to: 1.0,
            intensity: 0.0,
            integral: 0.0,
            frequency: 0.5,
            n_samples: 1,
        };
        assert!(f.mz.is_nan(), "NaN mz stays NaN");
        assert!(f.rt.is_infinite() && f.rt > 0.0, "positive infinity survives");
        assert!(f.from.is_infinite() && f.from < 0.0, "negative infinity survives");
        assert_eq!(f.to, 1.0, "finite values are unchanged");
    }

    #[test]
    fn test_feature_keeps_exact_values() {
        use crate::utilities::find_features::Feature;
        let f = Feature {
            mz: f64::NAN,
            rt: f64::INFINITY,
            from: 0.0,
            to: 1.0,
            intensity: 0.0,
            integral: 0.0,
            n_points: 0,
            noise: 0.0,
        };
        assert!(f.mz.is_nan(), "NaN mz stays NaN");
        assert!(f.rt.is_infinite(), "infinity survives");
        assert_eq!(f.to, 1.0, "finite values are unchanged");
    }

    #[test]
    fn test_peak_keeps_exact_values() {
        use crate::utilities::structs::Peak;
        let p = Peak {
            from: f64::NAN,
            to: f64::INFINITY,
            rt: 1.5,
            integral: 0.0,
            intensity: 0.0,
            n_points: 3,
            noise: 0.0,
            r2: None,
        };
        assert!(p.from.is_nan(), "NaN from stays NaN");
        assert!(p.to.is_infinite(), "infinity survives");
        assert_eq!(p.rt, 1.5, "finite rt unchanged");
        assert!(p.r2.is_none(), "absent r2 stays absent");
    }

    #[test]
    fn test_peak_default_is_all_zero() {
        use crate::utilities::structs::Peak;
        let p = Peak::default();
        assert_eq!(p.from, 0.0);
        assert_eq!(p.to, 0.0);
        assert_eq!(p.rt, 0.0);
        assert_eq!(p.n_points, 0);
        assert_eq!(p.noise, 0.0);
        assert!(p.r2.is_none());
    }

    #[test]
    fn test_spectrum_summary_unknown_keeps_not_a_number() {
        let summary = SpectrumSummary::unknown();
        assert!(summary.rt_seconds.is_nan(), "rt_seconds stays NaN");
        assert!(summary.base_peak_mz.is_nan(), "base_peak_mz stays NaN");
        assert!(summary.selected_ion_mz.is_nan(), "selected_ion_mz stays NaN");
        assert!(summary.base_peak_int.is_nan(), "base_peak_int stays NaN");
        assert!(
            summary.total_ion_current.is_nan(),
            "total_ion_current stays NaN"
        );
    }

    #[test]
    fn test_rt_shifted_samples_report_apex_mz() {
        let peak_mz = 200.0;
        let shifts = [4.90, 4.95, 5.00, 5.05, 5.10];

        let samples: Vec<Vec<u8>> = shifts
            .iter()
            .map(|&peak_rt| {
                let mzml = build_mzml_single_peak(peak_mz, peak_rt, 10_000.0, 4.0, 6.0, 60);
                to_ion_bytes(&mzml)
            })
            .collect();

        let dir = write_ion_dir("rt_shifted_apex", &samples);
        let features = get_features(
            dir.to_str().unwrap(),
            FromTo { from: 3.0, to: 7.0 },
            feature_cfg_for_mz(peak_mz),
            alignment_cfg(),
            1,
        )
        .expect("get_features should succeed");
        let _ = fs::remove_dir_all(&dir);

        let best = features.iter().max_by_key(|f| f.n_samples).unwrap();
        assert_eq!(best.n_samples, shifts.len(), "all shifted samples cluster");
        assert!(
            (best.mz - peak_mz).abs() < 0.001,
            "apex mz {:.5} should equal {peak_mz}",
            best.mz
        );
    }

    #[test]
    fn test_consensus_mz_is_deterministic_across_cores() {
        let peak_mz = 350.0;
        let shifts = [6.95, 7.00, 7.05, 6.90];

        let samples: Vec<Vec<u8>> = shifts
            .iter()
            .map(|&peak_rt| {
                to_ion_bytes(&build_mzml_single_peak(
                    peak_mz, peak_rt, 9_000.0, 6.0, 8.0, 60,
                ))
            })
            .collect();

        let run = |cores: usize| {
            let dir = write_ion_dir(&format!("determinism_{cores}"), &samples);
            let mut features = get_features(
                dir.to_str().unwrap(),
                FromTo { from: 5.0, to: 9.0 },
                feature_cfg_for_mz(peak_mz),
                alignment_cfg(),
                cores,
            )
            .expect("get_features should succeed");
            let _ = fs::remove_dir_all(&dir);
            features.sort_by(|a, b| a.mz.partial_cmp(&b.mz).unwrap());
            features.into_iter().map(|f| f.mz).collect::<Vec<f64>>()
        };

        let one = run(1);
        let many = run(4);
        assert_eq!(one.len(), many.len(), "same number of features");
        for (a, b) in one.iter().zip(many.iter()) {
            assert!(
                (a - b).abs() < 1e-9,
                "mz must be deterministic across cores: {a} vs {b}"
            );
        }
    }

    #[test]
    fn test_weighted_intensity_median_selection_reports_mz() {
        let peak_mz = 250.0;
        let peak_rt = 5.0;

        let samples: Vec<Vec<u8>> = (0..4)
            .map(|_| {
                to_ion_bytes(&build_mzml_single_peak(
                    peak_mz, peak_rt, 10_000.0, 4.0, 6.0, 60,
                ))
            })
            .collect();

        let dir = write_ion_dir("weighted_median_select", &samples);
        let mut cfg = alignment_cfg();
        cfg.mz_estimator = MzEstimatorKind::WeightedIntensityMedian;
        let features = get_features(
            dir.to_str().unwrap(),
            FromTo { from: 3.0, to: 7.0 },
            feature_cfg_for_mz(peak_mz),
            cfg,
            1,
        )
        .expect("get_features should succeed");
        let _ = fs::remove_dir_all(&dir);

        let best = features.iter().max_by_key(|f| f.n_samples).unwrap();
        assert!(
            (best.mz - peak_mz).abs() < 0.001,
            "weighted-median mz {:.5} should equal {peak_mz}",
            best.mz
        );
    }

    #[test]
    fn test_two_co_eluting_masses_split_into_two_features() {
        let mz_a = 200.0300;
        let mz_b = 200.0330;
        let peak_rt = 5.0;
        let n_samples = 6;

        let samples: Vec<Vec<u8>> = (0..n_samples)
            .map(|_| {
                to_ion_bytes(&build_mzml_two_peaks(
                    mz_a, mz_b, peak_rt, 10_000.0, 4.0, 6.0, 60,
                ))
            })
            .collect();

        let dir = write_ion_dir("co_eluting_split", &samples);
        let features = get_features(
            dir.to_str().unwrap(),
            FromTo { from: 3.0, to: 7.0 },
            feature_cfg_for_mz_range(mz_a, mz_b),
            alignment_cfg(),
            1,
        )
        .expect("get_features should succeed");
        let _ = fs::remove_dir_all(&dir);

        let reported: Vec<f64> = features.iter().map(|f| f.mz).collect();
        let near = |target: f64| {
            features
                .iter()
                .filter(|f| (f.mz - target).abs() < 0.001 && (f.rt - peak_rt).abs() < 0.2)
                .count()
        };
        assert_eq!(near(mz_a), 1, "mass A reported once, got {reported:?}");
        assert_eq!(near(mz_b), 1, "mass B reported once, got {reported:?}");
    }

    fn count_near(features: &[ConsensusFeature], mz: f64, rt: f64) -> usize {
        features
            .iter()
            .filter(|f| (f.mz - mz).abs() < 0.0008 && (f.rt - rt).abs() < 0.2)
            .count()
    }

    #[test]
    fn test_two_masses_within_cutoff_stay_one_feature() {
        let mz_a = 200.0300;
        let mz_b = 200.0310;
        let peak_rt = 5.0;
        let samples: Vec<Vec<u8>> = (0..6)
            .map(|_| {
                to_ion_bytes(&build_mzml_masses(
                    &[mz_a, mz_b],
                    peak_rt,
                    10_000.0,
                    4.0,
                    6.0,
                    60,
                ))
            })
            .collect();
        let dir = write_ion_dir("within_cutoff", &samples);
        let features = get_features(
            dir.to_str().unwrap(),
            FromTo { from: 3.0, to: 7.0 },
            feature_cfg_for_mz_range(mz_a, mz_b),
            alignment_cfg(),
            1,
        )
        .expect("get_features should succeed");
        let _ = fs::remove_dir_all(&dir);

        let near = features
            .iter()
            .filter(|f| f.mz > 200.028 && f.mz < 200.033 && (f.rt - peak_rt).abs() < 0.2)
            .count();
        assert_eq!(
            near,
            1,
            "masses closer than the cutoff must stay one feature, got {:?}",
            features.iter().map(|f| f.mz).collect::<Vec<_>>()
        );
    }

    // Regression: a single peak in PROFILE mode (one mass sampled at several
    // evenly-spaced m/z points, like TripleTOF data) must resolve to ONE mass,
    // not one per profile point. Currently `apex_scan_centroids` returns every
    // raw profile point, so resolve_cluster takes the multi-mass path and shatters
    // the peak. See mtbls733 m/z 325.177 @ rt 6.
    #[test]
    fn test_single_profile_peak_is_one_mass() {
        let center = 325.177;
        let step = 0.00253;
        let n_mz = 13;
        let mz_sigma = 0.005;
        let scans: Vec<(f64, Vec<f64>, Vec<f64>)> = (0..41)
            .map(|i| {
                let rt = 4.0 + 2.0 * i as f64 / 40.0;
                let rt_height = gaussian(rt, 5.0, 0.1, 100_000.0);
                let mut mzs = Vec::with_capacity(n_mz);
                let mut intensities = Vec::with_capacity(n_mz);
                for k in 0..n_mz {
                    let m = center + (k as f64 - (n_mz / 2) as f64) * step;
                    mzs.push(m);
                    intensities.push(rt_height * gaussian(m, center, mz_sigma, 1.0));
                }
                (rt, mzs, intensities)
            })
            .collect();
        let bounds = bounds_at(center, 5.0, 4.85, 5.15);
        let resolved = resolve_cluster(
            &scans,
            None,
            &bounds,
            center - 0.0065,
            center + 0.0065,
            20.0,
            None,
            SpectrumKind::Profile,
        );
        let masses = resolved.masses.expect("a real peak should measure a mass");
        let mzs: Vec<f64> = masses.iter().map(|m| m.mz).collect();
        assert_eq!(
            masses.len(),
            1,
            "one profile peak is a single mass, but it split into {} at m/z {:?}",
            masses.len(),
            mzs
        );
    }

    // Same defect end to end: one profile peak across several samples must be a
    // single consensus feature present in all of them, not a scatter of fragments.
    #[test]
    fn test_single_profile_peak_yields_one_feature() {
        let center = 325.177;
        let peak_rt = 5.0;
        let samples: Vec<Vec<u8>> = (0..6)
            .map(|_| {
                to_ion_bytes(&build_mzml_profile_peak(
                    center, 0.00253, 13, 0.005, peak_rt, 100_000.0, 4.0, 6.0, 60,
                ))
            })
            .collect();
        let dir = write_ion_dir("profile_single_mass", &samples);
        let features = get_features(
            dir.to_str().unwrap(),
            FromTo { from: 3.0, to: 7.0 },
            feature_cfg_for_mz(center),
            alignment_cfg(),
            1,
        )
        .expect("get_features should succeed");
        let _ = fs::remove_dir_all(&dir);

        let near: Vec<_> = features
            .iter()
            .filter(|f| (f.mz - center).abs() < 0.01 && (f.rt - peak_rt).abs() < 0.2)
            .collect();
        assert_eq!(
            near.len(),
            1,
            "one profile peak must be one feature, got {:?}",
            near.iter()
                .map(|f| (f.mz, f.rt, f.frequency))
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_three_co_eluting_masses_split_into_three() {
        let mzs = [200.030, 200.034, 200.038];
        let peak_rt = 5.0;
        let samples: Vec<Vec<u8>> = (0..6)
            .map(|_| to_ion_bytes(&build_mzml_masses(&mzs, peak_rt, 10_000.0, 4.0, 6.0, 60)))
            .collect();
        let dir = write_ion_dir("three_masses", &samples);
        let features = get_features(
            dir.to_str().unwrap(),
            FromTo { from: 3.0, to: 7.0 },
            feature_cfg_for_mz_range(mzs[0], mzs[2]),
            alignment_cfg(),
            1,
        )
        .expect("get_features should succeed");
        let _ = fs::remove_dir_all(&dir);

        let reported: Vec<f64> = features.iter().map(|f| f.mz).collect();
        for &mz in &mzs {
            assert_eq!(
                count_near(&features, mz, peak_rt),
                1,
                "mass {mz} once, got {reported:?}"
            );
        }
    }

    #[test]
    fn test_wider_tolerance_keeps_masses_merged() {
        let mz_a = 200.0300;
        let mz_b = 200.0330;
        let peak_rt = 5.0;
        let samples: Vec<Vec<u8>> = (0..6)
            .map(|_| {
                to_ion_bytes(&build_mzml_masses(
                    &[mz_a, mz_b],
                    peak_rt,
                    10_000.0,
                    4.0,
                    6.0,
                    60,
                ))
            })
            .collect();
        let dir = write_ion_dir("wide_tolerance", &samples);
        let mut cfg = alignment_cfg();
        cfg.eic_options.ppm_tolerance = 60.0;
        let features = get_features(
            dir.to_str().unwrap(),
            FromTo { from: 3.0, to: 7.0 },
            feature_cfg_for_mz_range(mz_a, mz_b),
            cfg,
            1,
        )
        .expect("get_features should succeed");
        let _ = fs::remove_dir_all(&dir);

        let near = features
            .iter()
            .filter(|f| f.mz > 200.028 && f.mz < 200.035 && (f.rt - peak_rt).abs() < 0.2)
            .count();
        assert_eq!(
            near,
            1,
            "a wider ppm tolerance must keep the masses merged, got {:?}",
            features.iter().map(|f| f.mz).collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_contaminant_without_peak_is_dropped() {
        let peak_mz = 200.0300;
        let flat_mz = 200.0340;
        let peak_rt = 5.0;
        let samples: Vec<Vec<u8>> = (0..6)
            .map(|_| {
                to_ion_bytes(&build_mzml_peak_and_flat(
                    peak_mz, flat_mz, peak_rt, 10_000.0, 4.0, 6.0, 60,
                ))
            })
            .collect();
        let dir = write_ion_dir("contaminant", &samples);
        let features = get_features(
            dir.to_str().unwrap(),
            FromTo { from: 3.0, to: 7.0 },
            feature_cfg_for_mz_range(peak_mz, flat_mz),
            alignment_cfg(),
            1,
        )
        .expect("get_features should succeed");
        let _ = fs::remove_dir_all(&dir);

        let reported: Vec<f64> = features.iter().map(|f| f.mz).collect();
        assert_eq!(
            count_near(&features, peak_mz, peak_rt),
            1,
            "the real mass is reported, got {reported:?}"
        );
        assert_eq!(
            count_near(&features, flat_mz, peak_rt),
            0,
            "a flat contaminant with no peak must be dropped, got {reported:?}"
        );
    }
    fn benchmark_clusterer() -> FeatureClusterer {
        FeatureClusterer {
            tolerance: MzTolerance {
                mz_absolute: 0.0025,
                ppm: 10.0,
            },
            rt_tolerance: 0.5,
        }
    }

    fn tag_rows(rows: &[(usize, f64, f64, f64, f64, f64)]) -> Vec<TaggedFeature> {
        rows.iter()
            .map(|(sample, mz, rt, from, to, intensity)| TaggedFeature {
                sample_index: *sample,
                feature: Feature {
                    mz: *mz,
                    rt: *rt,
                    from: *from,
                    to: *to,
                    intensity: *intensity,
                    ..Default::default()
                },
            })
            .collect()
    }

    fn count_features(clusters: &[Vec<TaggedFeature>]) -> usize {
        clusters.iter().map(|cluster| cluster.len()).sum()
    }

    fn rts_of(cluster: &[TaggedFeature]) -> Vec<f64> {
        let mut times: Vec<f64> = cluster.iter().map(|t| t.feature.rt).collect();
        times.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        times
    }

    fn rt_span_of(cluster: &[TaggedFeature]) -> f64 {
        let times = rts_of(cluster);
        times[times.len() - 1] - times[0]
    }

    fn mz_span_of(cluster: &[TaggedFeature]) -> f64 {
        let mut masses: Vec<f64> = cluster.iter().map(|t| t.feature.mz).collect();
        masses.sort_by(|a, b| a.partial_cmp(b).unwrap_or(Ordering::Equal));
        masses[masses.len() - 1] - masses[0]
    }

    fn count_samples_in(cluster: &[TaggedFeature]) -> usize {
        let mut samples: Vec<usize> = cluster.iter().map(|t| t.sample_index).collect();
        samples.sort_unstable();
        samples.dedup();
        samples.len()
    }

    fn find_cluster_holding_rt(
        clusters: &[Vec<TaggedFeature>],
        rt: f64,
    ) -> Option<&Vec<TaggedFeature>> {
        clusters
            .iter()
            .find(|cluster| cluster.iter().any(|t| (t.feature.rt - rt).abs() < 1e-9))
    }

    fn find_cluster_holding_mz(
        clusters: &[Vec<TaggedFeature>],
        mz: f64,
    ) -> Option<&Vec<TaggedFeature>> {
        clusters
            .iter()
            .find(|cluster| cluster.iter().any(|t| (t.feature.mz - mz).abs() < 1e-9))
    }

    fn partition_of(clusters: &[Vec<TaggedFeature>]) -> Vec<Vec<f64>> {
        let mut groups: Vec<Vec<f64>> = clusters.iter().map(|c| rts_of(c)).collect();
        groups.sort_by(|a, b| a[0].partial_cmp(&b[0]).unwrap_or(Ordering::Equal));
        groups
    }

    fn copy_of(cluster: &[TaggedFeature]) -> Vec<TaggedFeature> {
        cluster
            .iter()
            .map(|t| TaggedFeature {
                sample_index: t.sample_index,
                feature: t.feature.clone(),
            })
            .collect()
    }

    fn slots_of(cluster: &[TaggedFeature], n_samples: usize) -> Vec<Option<Feature>> {
        assign_best_per_sample(copy_of(cluster), n_samples)
    }

    fn count_samples_near(slots: &[Option<Feature>], rt: f64, allowed: f64) -> usize {
        slots
            .iter()
            .flatten()
            .filter(|f| (f.rt - rt).abs() <= allowed)
            .count()
    }

    fn count_features_near(features: &[TaggedFeature], rt: f64, allowed: f64) -> usize {
        features
            .iter()
            .filter(|t| (t.feature.rt - rt).abs() <= allowed)
            .count()
    }

    fn count_clusters_near(clusters: &[Vec<TaggedFeature>], rt: f64, allowed: f64) -> usize {
        clusters
            .iter()
            .filter(|cluster| {
                let mut times = rts_of(cluster);
                (median(&mut times) - rt).abs() <= allowed
            })
            .count()
    }

    fn rts_of_the_strongest(features: &[TaggedFeature], how_many: usize) -> Vec<f64> {
        let mut ranked: Vec<&TaggedFeature> = features.iter().collect();
        ranked.sort_by(|a, b| {
            b.feature
                .intensity
                .partial_cmp(&a.feature.intensity)
                .unwrap_or(Ordering::Equal)
        });
        ranked
            .into_iter()
            .take(how_many)
            .map(|t| t.feature.rt)
            .collect()
    }

    fn find_pairs_of_one_sample_that_do_not_overlap(cluster: &[TaggedFeature]) -> Vec<(f64, f64)> {
        let mut apart = Vec::new();
        for (position, one) in cluster.iter().enumerate() {
            for other in cluster.iter().skip(position + 1) {
                if one.sample_index != other.sample_index {
                    continue;
                }
                let shared =
                    one.feature.to.min(other.feature.to) - one.feature.from.max(other.feature.from);
                if shared <= 0.0 {
                    apart.push((one.feature.rt, other.feature.rt));
                }
            }
        }
        apart
    }

    fn get_features_at_278_12() -> Vec<TaggedFeature> {
        tag_rows(&[
            (2, 278.1247092, 1.507316, 1.306651, 1.623667, 6902308.1),
            (1, 278.1247344, 1.376137, 1.268542, 1.384857, 5148468.9),
            (0, 278.1247368, 1.167972, 1.136020, 1.240611, 4200416.7),
            (1, 278.1247559, 0.585820, 0.556843, 0.626382, 24165290.4),
            (2, 278.1247881, 0.580362, 0.557179, 0.629618, 19620234.3),
            (0, 278.1247936, 0.581682, 0.558518, 0.639630, 22771507.0),
            (1, 278.1248127, 1.222046, 1.073778, 1.268542, 3468844.7),
            (1, 278.1248132, 1.492446, 1.384857, 1.620422, 7335516.2),
            (0, 278.1248151, 1.275467, 1.240611, 1.284179, 5411123.4),
            (3, 278.1248154, 0.584686, 0.558601, 0.633942, 22301049.3),
            (4, 278.1248791, 1.497688, 1.328979, 1.602412, 6274111.2),
            (5, 278.1248933, 1.246061, 1.196641, 1.298394, 3275898.9),
            (0, 278.1249449, 1.377168, 1.284179, 1.388803, 6569379.3),
            (1, 278.1249474, 0.649577, 0.626382, 0.745360, 4136895.2),
            (7, 278.1249614, 1.250197, 1.133894, 1.302536, 3122101.7),
            (2, 278.1249713, 0.664416, 0.629618, 0.763167, 3952980.1),
            (9, 278.1249774, 1.207801, 0.888013, 1.228157, 2760609.1),
            (0, 278.1250003, 1.449824, 1.388803, 1.603906, 8389086.3),
            (2, 278.1250052, 1.277586, 1.204908, 1.283401, 3885649.2),
            (1, 278.1250113, 1.059252, 1.006937, 1.073778, 2935155.8),
            (1, 278.1250159, 1.001129, 0.893578, 1.006937, 2668988.8),
            (2, 278.1250182, 3.324813, 3.304351, 3.362783, 878420.6),
            (4, 278.1250226, 0.587815, 0.555941, 0.628407, 20144893.4),
            (2, 278.1250240, 0.946237, 0.882319, 0.972391, 2040532.4),
            (0, 278.1250249, 0.668613, 0.639630, 0.761458, 4300743.0),
            (8, 278.1250308, 0.584550, 0.561376, 0.639619, 22074602.5),
            (3, 278.1250553, 0.744277, 0.735566, 0.758802, 909908.0),
            (4, 278.1250623, 1.439497, 1.328979, 1.602412, 6558257.0),
            (6, 278.1250666, 1.558022, 1.552203, 1.648223, 5121739.6),
            (5, 278.1250763, 0.583666, 0.557587, 0.635845, 21310584.0),
            (3, 278.1250772, 0.665838, 0.633942, 0.721046, 4606471.2),
            (3, 278.1250783, 1.122097, 0.968017, 1.142445, 2825406.4),
            (4, 278.1250863, 0.671907, 0.628407, 0.741614, 3669121.5),
            (3, 278.1250891, 0.965115, 0.895373, 0.968017, 1588326.2),
            (4, 278.1250955, 1.279525, 1.186455, 1.328979, 4676465.0),
            (9, 278.1250980, 1.533613, 1.399766, 1.542336, 7161487.5),
            (3, 278.1250981, 1.200581, 1.142445, 1.235488, 3609348.4),
            (0, 278.1251151, 0.947261, 0.894988, 0.961771, 2412317.5),
            (0, 278.1251182, 1.005325, 0.961771, 1.080812, 3118775.4),
            (7, 278.1251202, 0.584918, 0.564648, 0.634191, 23356206.7),
            (6, 278.1251275, 1.275876, 1.165368, 1.319513, 3663239.5),
            (4, 278.1251291, 1.160286, 1.105032, 1.186455, 3467646.9),
            (9, 278.1251355, 0.577370, 0.559982, 0.644021, 22202370.8),
            (4, 278.1251452, 3.333497, 3.304262, 3.359797, 821505.1),
            (4, 278.1251470, 1.038156, 0.988708, 1.046879, 2644703.4),
            (5, 278.1251506, 1.528174, 1.516545, 1.638769, 4866954.3),
            (3, 278.1251554, 1.407055, 1.340178, 1.415774, 4971958.6),
            (6, 278.1251634, 0.587283, 0.561214, 0.630747, 22864950.4),
            (5, 278.1252031, 0.653243, 0.635845, 0.728711, 3409901.9),
            (5, 278.1252163, 1.359465, 1.196641, 1.365283, 4159205.9),
            (9, 278.1252195, 1.274687, 1.129295, 1.283410, 3077606.2),
            (9, 278.1252491, 1.370676, 1.283410, 1.399766, 3575587.3),
            (7, 278.1252505, 0.654487, 0.634191, 0.756125, 2883346.7),
            (5, 278.1252524, 1.068746, 0.882751, 1.074567, 2564009.8),
            (5, 278.1252613, 3.343977, 3.291307, 3.364473, 975462.9),
            (5, 278.1252674, 1.487451, 1.478724, 1.638769, 6779748.6),
            (6, 278.1252901, 0.778773, 0.770063, 0.787484, 610950.7),
            (8, 278.1253026, 1.174429, 1.078484, 1.203512, 2901169.5),
            (8, 278.1253057, 0.665725, 0.639619, 0.747091, 2865672.8),
            (6, 278.1253068, 0.694584, 0.630747, 0.711999, 3116408.9),
            (3, 278.1253081, 1.328545, 1.305283, 1.340178, 4071073.0),
            (7, 278.1253313, 3.333886, 3.307546, 3.357281, 1266488.9),
            (7, 278.1253349, 0.936270, 0.886877, 0.939178, 1927652.6),
            (6, 278.1253419, 0.999674, 0.982232, 1.078161, 2355544.6),
            (7, 278.1253450, 1.130985, 0.939178, 1.133894, 2264437.1),
            (9, 278.1253452, 0.946142, 0.888013, 0.981031, 1875775.5),
            (6, 278.1253469, 1.142114, 1.078161, 1.165368, 3008511.6),
            (9, 278.1253475, 1.013008, 0.981031, 1.062425, 2043276.6),
            (9, 278.1253571, 1.103128, 1.085681, 1.129295, 2673913.8),
            (8, 278.1253683, 3.333824, 3.295795, 3.351367, 908105.5),
            (6, 278.1253686, 1.499828, 1.319513, 1.648223, 6592798.8),
            (9, 278.1253704, 1.667491, 1.658756, 1.702420, 695893.2),
            (8, 278.1253742, 0.999999, 0.889529, 1.031972, 1748936.7),
            (8, 278.1253848, 1.072672, 1.031972, 1.078484, 2310676.5),
            (5, 278.1254003, 1.141416, 1.074567, 1.196641, 3544183.9),
            (2, 278.1254398, 1.216541, 1.132233, 1.231078, 3874205.3),
            (8, 278.1254562, 1.410007, 1.325641, 1.412917, 4728272.1),
            (7, 278.1254721, 1.532324, 1.302536, 1.674923, 6413670.3),
            (8, 278.1254853, 1.517662, 1.412917, 1.671884, 6729067.0),
            (8, 278.1254961, 1.305287, 1.203512, 1.325641, 3387420.0),
            (9, 278.1255480, 0.652732, 0.644021, 0.751423, 3984235.6),
            (9, 278.1255900, 1.481226, 1.228157, 1.487044, 5431401.8),
            (4, 278.1263846, 0.918951, 0.875338, 0.924761, 1442638.4),
        ])
    }

    fn get_two_masses_at_304_20() -> Vec<TaggedFeature> {
        tag_rows(&[
            (3, 304.2025617, 3.288708, 2.073540, 3.376424, 49974342.6),
            (1, 304.2025916, 3.286318, 2.196820, 3.400336, 58018602.7),
            (2, 304.2025995, 3.286782, 2.080718, 3.438820, 49202635.2),
            (0, 304.2026136, 3.285310, 2.281811, 3.369970, 55261380.6),
            (4, 304.2026992, 3.280833, 2.289790, 3.371510, 50700216.7),
            (6, 304.2027710, 3.283394, 2.332799, 3.371198, 47752762.0),
            (8, 304.2028932, 3.278183, 2.164176, 3.453827, 43592643.8),
            (9, 304.2029127, 3.279758, 2.174331, 3.411447, 43457626.2),
            (5, 304.2029461, 3.282526, 2.195044, 3.455203, 48777155.5),
            (7, 304.2029671, 3.281172, 2.222525, 3.439274, 53475197.0),
            (3, 304.2055435, 3.514121, 3.464277, 4.582163, 66571799.6),
            (5, 304.2056400, 3.516797, 3.472782, 3.634202, 70860810.5),
            (0, 304.2056454, 3.501548, 3.466423, 3.636235, 75983589.0),
            (4, 304.2056560, 3.509285, 3.459392, 4.476941, 68828069.1),
            (8, 304.2057315, 3.533049, 3.483158, 4.632856, 66894172.2),
        ])
    }

    fn get_two_peaks_at_334_06() -> Vec<TaggedFeature> {
        tag_rows(&[
            (1, 334.0648000, 4.978775, 4.904848, 5.201779, 32077699.3),
            (1, 334.0649065, 4.422416, 4.357508, 4.552294, 18369310.0),
            (2, 334.0649124, 4.426315, 4.361366, 4.500165, 14973467.7),
            (2, 334.0649196, 4.983409, 4.921268, 5.093314, 28952532.4),
            (0, 334.0649864, 4.943248, 4.872209, 5.147846, 28192380.7),
            (0, 334.0650187, 4.393444, 4.334567, 4.464108, 18577135.8),
            (3, 334.0651103, 4.425401, 4.372260, 4.511034, 15989103.1),
            (3, 334.0651231, 4.994339, 4.929229, 5.110301, 28946592.1),
            (4, 334.0651420, 5.004766, 4.936712, 5.237302, 25614467.4),
            (5, 334.0652053, 4.482000, 4.428681, 4.574185, 17172780.7),
            (4, 334.0652747, 4.432620, 4.376525, 4.518378, 16894537.9),
            (8, 334.0652997, 4.576442, 4.519730, 4.665418, 18369066.8),
            (5, 334.0653824, 4.591966, 4.580105, 4.600872, 251859.3),
            (6, 334.0653903, 4.515814, 4.462254, 4.602049, 17279632.9),
            (9, 334.0654215, 4.613515, 4.550994, 4.705320, 17045829.1),
            (7, 334.0654704, 4.549616, 4.492919, 4.635714, 16668725.4),
        ])
    }

    fn get_one_compound_at_334_06() -> Vec<TaggedFeature> {
        tag_rows(&[
            (1, 334.0649065, 4.422416, 4.357508, 4.552294, 18369310.0),
            (2, 334.0649124, 4.426315, 4.361366, 4.500165, 14973467.7),
            (0, 334.0650187, 4.393444, 4.334567, 4.464108, 18577135.8),
            (3, 334.0651103, 4.425401, 4.372260, 4.511034, 15989103.1),
            (5, 334.0652053, 4.482000, 4.428681, 4.574185, 17172780.7),
            (4, 334.0652747, 4.432620, 4.376525, 4.518378, 16894537.9),
            (8, 334.0652997, 4.576442, 4.519730, 4.665418, 18369066.8),
            (5, 334.0653824, 4.591966, 4.580105, 4.600872, 251859.3),
            (6, 334.0653903, 4.515814, 4.462254, 4.602049, 17279632.9),
            (9, 334.0654215, 4.613515, 4.550994, 4.705320, 17045829.1),
            (7, 334.0654704, 4.549616, 4.492919, 4.635714, 16668725.4),
        ])
    }

    fn get_half_the_samples_at_334_06() -> Vec<TaggedFeature> {
        tag_rows(&[
            (5, 334.0652296, 5.082625, 5.008229, 5.312290, 26210077.0),
            (6, 334.0652362, 5.116502, 5.045162, 5.229948, 27320614.4),
            (8, 334.0653165, 5.203912, 5.135477, 5.311684, 27176893.4),
            (7, 334.0654226, 5.161793, 5.084573, 5.272240, 26995653.8),
            (9, 334.0654883, 5.252645, 5.172190, 5.456673, 28589492.1),
        ])
    }

    fn get_a_chain_at_344_85() -> Vec<TaggedFeature> {
        tag_rows(&[
            (2, 344.8501502, 35.954479, 35.945719, 35.963250, 150645.9),
            (1, 344.8501695, 3.482461, 3.473619, 3.491271, 108626.1),
            (8, 344.8505195, 31.456582, 31.447419, 31.465751, 40287.9),
            (5, 344.8505628, 2.536373, 2.527621, 2.545125, 233944.4),
            (5, 344.8520941, 32.999978, 32.990576, 33.009947, 35838.7),
            (5, 344.8521588, 33.533814, 33.524528, 33.543122, 47844.9),
            (7, 344.8521848, 33.274560, 33.265180, 33.283945, 35044.5),
            (7, 344.8522092, 32.934916, 32.925545, 32.944341, 36072.1),
            (7, 344.8522210, 34.249136, 34.239625, 34.258658, 39052.7),
            (9, 344.8522290, 33.336877, 33.327515, 33.346198, 56316.1),
            (7, 344.8538670, 33.903642, 33.894147, 33.916332, 42367.8),
            (8, 344.8538746, 33.633224, 33.623876, 33.642528, 62789.2),
            (5, 344.8538820, 35.145525, 35.136283, 35.154810, 58640.3),
            (9, 344.8538838, 33.591852, 33.581896, 33.601227, 54028.3),
            (8, 344.8539099, 35.416671, 35.407205, 35.426155, 35921.8),
            (8, 344.8539148, 33.144597, 33.135174, 33.154034, 40619.5),
            (8, 344.8539156, 34.110569, 34.101083, 34.120087, 33945.4),
            (8, 344.8539399, 34.631640, 34.622146, 34.641129, 36282.8),
            (8, 344.8539651, 34.812749, 34.803249, 34.822230, 36173.0),
            (7, 344.8539776, 22.949721, 22.939480, 22.959842, 22663.0),
            (2, 344.8552176, 0.524357, 0.515444, 0.533296, 84139.2),
            (1, 344.8552891, 22.087077, 22.077577, 22.097185, 48163.4),
            (3, 344.8554674, 3.042745, 3.033966, 3.051524, 181699.5),
            (5, 344.8555812, 33.097748, 33.088334, 33.107144, 39023.0),
            (9, 344.8556292, 33.516467, 33.507092, 33.525800, 48482.7),
            (7, 344.8556444, 33.723314, 33.713919, 33.732752, 47847.3),
            (9, 344.8557172, 32.516885, 32.507935, 32.525908, 85546.3),
            (1, 344.8569989, 32.058797, 32.049814, 32.067826, 92887.2),
        ])
    }

    fn get_one_lane_at_334_07() -> Vec<TaggedFeature> {
        tag_rows(&[
            (2, 334.0698058, 0.571683, 0.562994, 0.580362, 927328.1),
            (3, 334.0699029, 0.578899, 0.570213, 0.587580, 967090.1),
            (1, 334.0716456, 6.866517, 6.736422, 7.320749, 21631986.3),
            (2, 334.0717092, 6.864335, 6.736871, 7.342225, 21994707.6),
            (5, 334.0717756, 5.686416, 5.611223, 5.759286, 406934.9),
            (0, 334.0718069, 6.886967, 6.766045, 7.690180, 29184679.5),
            (4, 334.0718946, 6.851507, 6.830094, 6.869783, 24712983.9),
            (6, 334.0719088, 5.703120, 5.654656, 5.800163, 495376.9),
            (5, 334.0719728, 6.818460, 6.713761, 7.777675, 18856959.9),
            (9, 334.0720983, 6.800715, 6.669114, 7.647857, 20902662.8),
            (6, 334.0721184, 6.808458, 6.691796, 6.829412, 18171095.1),
            (8, 334.0721202, 6.815308, 6.686901, 6.979176, 19082872.5),
            (7, 334.0733961, 5.761151, 5.688431, 5.830771, 1377825.7),
            (8, 334.0734113, 5.801595, 5.716654, 5.886305, 1407823.8),
            (9, 334.0736823, 5.832974, 5.763242, 5.920607, 1309727.5),
            (0, 334.0750367, 6.883946, 6.766045, 7.760450, 216373810.0),
            (9, 334.0751879, 6.809771, 6.666114, 7.656988, 152095805.5),
            (5, 334.0752499, 6.818460, 6.710758, 7.825896, 136867551.5),
            (6, 334.0752562, 6.811445, 6.688792, 7.027477, 121333751.8),
            (7, 334.0753027, 6.811825, 6.689390, 7.706503, 129595366.0),
            (8, 334.0753878, 6.809296, 6.683895, 7.657127, 142263580.5),
            (1, 334.0764831, 6.872575, 6.736422, 7.025099, 203144611.5),
            (2, 334.0765484, 6.849118, 6.733859, 7.811690, 207865523.8),
            (3, 334.0766303, 6.879052, 6.745378, 7.836011, 200898625.2),
            (4, 334.0767804, 6.845409, 6.726789, 6.875886, 193744980.0),
        ])
    }

    fn get_one_lane_at_140_01() -> Vec<TaggedFeature> {
        tag_rows(&[
            (1, 140.0107840, 29.942248, 3.837268, 33.191101, 8171617.9),
            (0, 140.0108165, 29.946965, 3.431270, 33.323329, 10102285.3),
            (2, 140.0108647, 29.934197, 3.421246, 33.253766, 8947355.5),
            (3, 140.0108725, 29.941619, 4.004275, 32.476754, 8846932.0),
            (4, 140.0108979, 29.920938, 3.714570, 31.915842, 8584938.0),
            (5, 140.0109527, 5.698514, 5.154134, 6.365001, 4533671.4),
            (6, 140.0109640, 29.923522, 21.371750, 31.773964, 8005937.3),
            (8, 140.0109658, 5.659150, 4.763549, 18.938066, 3776765.0),
            (7, 140.0109690, 29.939499, 4.644584, 31.861472, 8011010.4),
            (5, 140.0109790, 27.437446, 23.039894, 28.972227, 4536703.3),
            (9, 140.0109885, 29.931927, 4.761795, 31.843057, 7215614.2),
            (5, 140.0109902, 29.928420, 28.972227, 31.724873, 6170569.3),
            (6, 140.0110016, 5.657673, 5.104616, 11.435238, 4073573.8),
            (8, 140.0110351, 29.941979, 18.938066, 31.800465, 7751646.8),
        ])
    }
    fn get_three_features_on_one_mass_ladder() -> Vec<TaggedFeature> {
        tag_rows(&[
            (0, 100.000, 1.0, 0.95, 1.05, 500.0),
            (1, 100.009, 5.0, 4.95, 5.05, 400.0),
            (2, 100.017, 1.0, 0.95, 1.05, 300.0),
        ])
    }

    fn get_identical_features_from_every_sample() -> Vec<TaggedFeature> {
        (0..10)
            .map(|sample| {
                tag_rows(&[(sample, 278.12512, 0.5847, 0.5401, 0.6293, 22000000.0)])
                    .pop()
                    .unwrap()
            })
            .collect()
    }

    #[test]
    fn grouping_keeps_every_feature_it_is_given() {
        let features = get_features_at_278_12();
        let given = features.len();
        let clusters = benchmark_clusterer().cluster(features);
        let kept = count_features(&clusters);
        assert_eq!(
            kept, given,
            "grouping must not delete features, {given} went in and {kept} came out"
        );
    }

    #[test]
    fn grouping_keeps_a_cluster_where_the_features_are() {
        let features = get_features_at_278_12();
        let detected = count_features_near(&features, 0.58, 0.3);
        let clusters = benchmark_clusterer().cluster(features);
        let found = count_clusters_near(&clusters, 0.58, 0.3);
        assert!(
            found > 0,
            "{detected} features sit at 0.58 min and grouping returned {found} clusters there"
        );
    }

    #[test]
    fn grouping_keeps_the_most_intense_features() {
        let features = get_features_at_278_12();
        let strongest = rts_of_the_strongest(&features, 10);
        let clusters = benchmark_clusterer().cluster(features);
        let lost: Vec<f64> = strongest
            .into_iter()
            .filter(|rt| find_cluster_holding_rt(&clusters, *rt).is_none())
            .collect();
        assert!(
            lost.is_empty(),
            "grouping deleted {} of the ten most intense features, at {:?}",
            lost.len(),
            lost
        );
    }

    #[test]
    fn three_features_are_enough_to_lose_one() {
        let features = get_three_features_on_one_mass_ladder();
        let given = features.len();
        let clusters = default_clusterer().cluster(features);
        let kept = count_features(&clusters);
        assert_eq!(
            kept,
            given,
            "{given} features were handed to grouping and {kept} came back, holding {:?}",
            clusters
                .iter()
                .map(|c| c.iter().map(|t| t.feature.mz).collect::<Vec<f64>>())
                .collect::<Vec<Vec<f64>>>()
        );
    }

    #[test]
    fn a_cluster_of_two_can_still_lose_a_member() {
        let clusters = default_clusterer().cluster(get_three_features_on_one_mass_ladder());
        let holding = find_cluster_holding_mz(&clusters, 100.017)
            .expect("the feature at 100.017 is in no cluster");
        assert_eq!(
            holding.len(),
            2,
            "100.000 and 100.017 were grouped together and the cluster came back holding {}",
            holding.len()
        );
    }

    #[test]
    fn a_chain_of_features_keeps_every_one_of_them() {
        let features = get_a_chain_at_344_85();
        let given = features.len();
        let clusters = benchmark_clusterer().cluster(features);
        let kept = count_features(&clusters);
        assert_eq!(
            kept, given,
            "a chain of {given} features came back holding {kept}"
        );
    }

    #[test]
    fn features_pushed_out_of_a_cluster_start_a_new_one() {
        let clusters = benchmark_clusterer().cluster(get_a_chain_at_344_85());
        assert!(
            find_cluster_holding_rt(&clusters, 32.058797).is_some(),
            "the feature at 32.0588 is in no cluster, grouping returned {:?}",
            partition_of(&clusters)
        );
    }

    #[test]
    fn the_stronger_mass_does_not_take_over_the_weaker_one() {
        let clusters = benchmark_clusterer().cluster(get_two_masses_at_304_20());
        let cluster = find_cluster_holding_rt(&clusters, 3.285310)
            .expect("the feature at 3.2853 is in no cluster");
        let held = count_samples_near(&slots_of(cluster, 10), 3.28, 0.05);
        assert_eq!(
            held,
            10,
            "304.2026 at 3.28 was found in ten samples and a cluster of {} features \
             leaves it in {held}",
            cluster.len()
        );
    }

    #[test]
    fn a_stronger_neighbour_does_not_evict_the_weaker_feature_of_the_same_sample() {
        let clusters = benchmark_clusterer().cluster(get_two_peaks_at_334_06());
        let cluster = find_cluster_holding_rt(&clusters, 4.432620)
            .expect("the feature at 4.4326 is in no cluster");
        let held = count_samples_near(&slots_of(cluster, 10), 4.5, 0.15);
        assert_eq!(
            held, 10,
            "334.0650 at 4.48 was found in ten samples and the stronger peak at 4.98 \
             leaves it in {held}"
        );
    }

    #[test]
    fn search_bounds_centre_on_the_compound_not_on_its_neighbour() {
        let clusters = benchmark_clusterer().cluster(get_two_peaks_at_334_06());
        let cluster = find_cluster_holding_rt(&clusters, 4.432620)
            .expect("the feature at 4.4326 is in no cluster");
        let bounds =
            compute_search_bounds(&slots_of(cluster, 10), 0.5).expect("the cluster filled no slot");
        assert!(
            (bounds.center_rt - 4.48).abs() <= 0.15,
            "the compound sits at 4.48 and the search centre landed at {}",
            bounds.center_rt
        );
    }

    #[test]
    fn a_cluster_never_holds_two_peaks_of_one_sample_that_do_not_overlap() {
        for cluster in benchmark_clusterer().cluster(get_two_peaks_at_334_06()) {
            let apart = find_pairs_of_one_sample_that_do_not_overlap(&cluster);
            assert!(
                apart.is_empty(),
                "one cluster holds {} pairs from a single sample whose peaks never overlap, \
                 the first is {:?}",
                apart.len(),
                apart[0]
            );
        }
    }

    #[test]
    fn a_compound_stays_in_one_cluster_when_a_weaker_ion_shares_its_time() {
        let clusters = benchmark_clusterer().cluster(get_one_lane_at_334_07());
        let cluster = find_cluster_holding_mz(&clusters, 334.0750367)
            .expect("the feature at 334.07504 is in no cluster");
        let held = cluster.iter().filter(|t| t.feature.mz >= 334.0740).count();
        assert_eq!(
            held,
            10,
            "334.0750 to 334.0768 is one compound in ten samples spanning {} in mass, \
             grouping left {held} of its features in this cluster",
            mz_span_of(cluster)
        );
    }

    #[test]
    fn a_weaker_ion_stays_out_of_the_compound_that_shares_its_time() {
        let clusters = benchmark_clusterer().cluster(get_one_lane_at_334_07());
        let cluster = find_cluster_holding_mz(&clusters, 334.0750367)
            .expect("the feature at 334.07504 is in no cluster");
        let intruders = cluster.iter().filter(|t| t.feature.mz < 334.0740).count();
        assert_eq!(
            intruders, 0,
            "334.0719 and 334.0750 are two ions at one retention time and the cluster \
             of the stronger one swallowed {intruders} features of the weaker"
        );
    }

    #[test]
    fn one_compound_with_a_wide_gap_stays_in_one_cluster() {
        let features = get_one_compound_at_334_06();
        let given = features.len();
        let clusters = benchmark_clusterer().cluster(features);
        assert_eq!(
            count_features(&clusters),
            given,
            "grouping must keep every feature it was given, returned {:?}",
            partition_of(&clusters)
        );
        let biggest = clusters.iter().map(|cluster| cluster.len()).max().unwrap();
        assert_eq!(
            biggest,
            10,
            "the ten real peaks, one per sample, belong in one cluster and grouping returned {:?}",
            partition_of(&clusters)
        );
    }

    #[test]
    fn a_compound_seen_in_half_the_samples_survives() {
        let features = get_half_the_samples_at_334_06();
        let given = features.len();
        let clusters = benchmark_clusterer().cluster(features);
        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0].len(), given);
        assert_eq!(count_samples_in(&clusters[0]), 5);
    }

    #[test]
    fn identical_features_from_different_samples_form_one_cluster() {
        let clusters = benchmark_clusterer().cluster(get_identical_features_from_every_sample());
        assert_eq!(clusters.len(), 1);
        assert_eq!(clusters[0].len(), 10);
    }

    #[test]
    fn a_lane_spread_over_the_whole_run_never_becomes_one_cluster() {
        let clusterer = benchmark_clusterer();
        let allowed = 2.0 * clusterer.rt_tolerance;
        for cluster in clusterer.cluster(get_one_lane_at_140_01()) {
            let span = rt_span_of(&cluster);
            assert!(
                span <= allowed,
                "one cluster of {} features from {} samples spans {span} min \
                 where {allowed} is allowed",
                cluster.len(),
                count_samples_in(&cluster)
            );
        }
    }

    fn read_per_sample_features(folder: &std::path::Path) -> Vec<TaggedFeature> {
        let mut names: Vec<std::path::PathBuf> = fs::read_dir(folder)
            .expect("no per-sample folder")
            .filter_map(|entry| entry.ok().map(|e| e.path()))
            .filter(|path| path.extension().and_then(|s| s.to_str()) == Some("tsv"))
            .collect();
        names.sort();

        let mut tagged = Vec::new();
        for (sample_index, path) in names.iter().enumerate() {
            let text = fs::read_to_string(path).expect("unreadable table");
            for line in text.lines().skip(1) {
                let cells: Vec<&str> = line.split('\t').collect();
                if cells.len() < 7 {
                    continue;
                }
                let value = |index: usize| cells[index].parse::<f64>().unwrap_or(f64::NAN);
                tagged.push(TaggedFeature {
                    sample_index,
                    feature: Feature {
                        mz: value(0),
                        rt: value(1),
                        from: value(2),
                        to: value(3),
                        intensity: value(4),
                        integral: value(5),
                        n_points: cells[6].parse::<usize>().unwrap_or(0),
                        noise: 0.0,
                    },
                });
            }
        }
        tagged
    }

    fn read_truth_table(path: &std::path::Path) -> Vec<(String, f64, f64)> {
        let text = fs::read_to_string(path).expect("no truth table");
        let mut rows = Vec::new();
        for line in text.lines().skip(1) {
            let cells: Vec<&str> = line.split('\t').collect();
            if cells.len() < 3 {
                continue;
            }
            match (cells[1].parse::<f64>(), cells[2].parse::<f64>()) {
                (Ok(mz), Ok(rt)) => rows.push((cells[0].to_string(), mz, rt)),
                _ => continue,
            }
        }
        rows
    }

    #[test]
    #[ignore = "needs debug/find_features, run with: cargo test --lib -- --ignored --nocapture"]
    fn grouping_recovers_the_truth_over_the_real_samples() {
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("no repository root")
            .join("debug");
        let tagged = read_per_sample_features(&root.join("find_features"));
        let truth = read_truth_table(&root.join("truth-qehc.tsv"));
        let samples = tagged
            .iter()
            .map(|item| item.sample_index)
            .max()
            .map(|last| last + 1)
            .unwrap_or(0);
        let given = tagged.len();
        println!("{given} features from {samples} samples, {} truth compounds", truth.len());

        let began = std::time::Instant::now();
        let clusters = benchmark_clusterer().cluster(tagged);
        let took = began.elapsed();
        let kept: usize = clusters.iter().map(|cluster| cluster.len()).sum();

        let mut summary: Vec<(f64, f64, usize)> = clusters
            .iter()
            .map(|cluster| {
                let mut masses: Vec<f64> = cluster.iter().map(|t| t.feature.mz).collect();
                let mut times: Vec<f64> = cluster.iter().map(|t| t.feature.rt).collect();
                let mut seen: Vec<usize> = cluster.iter().map(|t| t.sample_index).collect();
                seen.sort_unstable();
                seen.dedup();
                (median(&mut masses), median(&mut times), seen.len())
            })
            .collect();
        summary.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));

        let masses: Vec<f64> = summary.iter().map(|row| row.0).collect();
        let mut complete = 0;
        let mut anywhere = 0;
        for (_, mz, rt) in &truth {
            let half = mz * 10.0 * 1e-6;
            let low = masses.partition_point(|value| *value < mz - half);
            let high = masses.partition_point(|value| *value <= mz + half);
            let near: Vec<&(f64, f64, usize)> = summary[low..high]
                .iter()
                .filter(|row| (row.1 - rt).abs() <= 0.3)
                .collect();
            if !near.is_empty() {
                anywhere += 1;
            }
            if near.iter().any(|row| row.2 == samples) {
                complete += 1;
            }
        }

        println!("clusters              : {}", clusters.len());
        println!("features in            : {given}");
        println!("features out           : {kept}");
        println!("features lost          : {}", given - kept);
        println!("truth found anywhere   : {anywhere} / {}", truth.len());
        println!("truth in all {samples} samples : {complete} / {}", truth.len());
        println!("elapsed                : {:.1} s", took.as_secs_f64());

        assert_eq!(given, kept, "grouping lost features");
    }
}
