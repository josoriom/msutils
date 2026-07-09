#[cfg(test)]
mod tests {
    use std::fs;
    use std::sync::Arc;

    use ionic::{ion::encode, mzml::structs::*};

    use crate::utilities::{
        calculate_eic::{CentroidScan, EicOptions, SpectrumKind, SpectrumSummary},
        find_features::{Feature, FindFeaturesOptions, MzScanGrid, MzTolerance},
        get_features::{
            AlignmentOptions, ConsensusFeature, FeatureClusterer, MzRtCluster, SearchBounds,
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
    fn test_growing_cluster_single_item() {
        let gc = MzRtCluster::new(make_tagged(0, 100.0, 1.0, 500.0));
        assert_eq!(gc.cached_median_mz, 100.0);
        assert_eq!(gc.cached_median_rt, 1.0);
    }

    #[test]
    fn test_growing_cluster_median_updates_on_push() {
        let mut gc = MzRtCluster::new(make_tagged(0, 100.0, 1.0, 500.0));
        gc.push(make_tagged(1, 102.0, 2.0, 300.0));
        assert_eq!(gc.cached_median_mz, 102.0);
    }

    #[test]
    fn test_growing_cluster_three_items_median() {
        let mut gc = MzRtCluster::new(make_tagged(0, 100.0, 1.0, 500.0));
        gc.push(make_tagged(1, 102.0, 2.0, 300.0));
        gc.push(make_tagged(2, 104.0, 3.0, 200.0));
        assert_eq!(gc.cached_median_mz, 102.0);
    }

    #[test]
    fn test_growing_cluster_four_items_upper_index() {
        let mut gc = MzRtCluster::new(make_tagged(0, 100.0, 1.0, 500.0));
        gc.push(make_tagged(1, 102.0, 2.0, 300.0));
        gc.push(make_tagged(2, 104.0, 3.0, 200.0));
        gc.push(make_tagged(3, 106.0, 4.0, 100.0));
        assert_eq!(gc.cached_median_mz, 104.0);
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
            make_tagged(0, 100.0, 1.0, 500.0),
            make_tagged(1, 100.005, 1.05, 400.0),
            make_tagged(0, 200.0, 2.0, 300.0),
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
            make_tagged(0, 100.0, 1.0, 500.0),
            make_tagged(1, 100.005, 1.0, 400.0),
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
            make_tagged(0, 100.0, 1.0, 500.0),
            make_tagged(1, 100.005, 1.05, 400.0),
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
        encode(mzml, 0, false, &mut bytes).expect("ion encode should succeed");
        bytes
    }

    fn write_ion_dir(tag: &str, samples: &[Vec<u8>]) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!("msutils_test_{tag}"));
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
        let resolved =
            resolve_cluster(&scans, None, &bounds, 99.99, 100.01, 20.0, None, SpectrumKind::Centroid);
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
        let resolved =
            resolve_cluster(&scans, None, &bounds, 299.99, 300.01, 20.0, None, SpectrumKind::Centroid);

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
            "filled feature mz should be set from the apex"
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
    fn test_nan_f64_serialises_as_zero_in_consensus_feature() {
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
        let json = serde_json::to_string(&f).expect("serialise ok");
        // Non-finite mz, rt, from must become 0, not null.
        assert!(
            !json.contains("null"),
            "non-finite fields must not produce null: {json}"
        );
        assert!(json.contains("\"mz\":0"), "NaN mz → 0: {json}");
        assert!(json.contains("\"rt\":0"), "+Inf rt → 0: {json}");
        assert!(json.contains("\"from\":0"), "-Inf from → 0: {json}");
        assert!(json.contains("\"to\":1.0"), "finite to unchanged: {json}");
    }

    #[test]
    fn test_nan_f64_serialises_as_zero_in_feature() {
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
        let json = serde_json::to_string(&f).expect("serialise ok");
        assert!(
            !json.contains("null"),
            "non-finite fields must not produce null: {json}"
        );
        assert!(json.contains("\"mz\":0"), "NaN mz → 0: {json}");
        assert!(json.contains("\"rt\":0"), "+Inf rt → 0: {json}");
    }

    #[test]
    fn test_nan_f64_serialises_as_zero_in_peak() {
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
        let json = serde_json::to_string(&p).expect("serialise ok");
        assert!(
            !json.contains("null"),
            "non-finite fields must not produce null: {json}"
        );
        assert!(json.contains("\"from\":0"), "NaN from → 0: {json}");
        assert!(json.contains("\"to\":0"), "+Inf to → 0: {json}");
        assert!(json.contains("\"rt\":1.5"), "finite rt unchanged: {json}");
        assert!(!json.contains("\"r2\""), "r2 must be skipped: {json}");
    }

    #[test]
    fn test_peak_default_serialises_as_zero_object() {
        use crate::utilities::structs::Peak;
        let json = serde_json::to_string(&Peak::default()).expect("serialise ok");
        assert!(json.starts_with('{'), "must be an object, not null: {json}");
        assert!(json.contains("\"from\":0"), "from defaults to 0: {json}");
        assert!(json.contains("\"to\":0"), "to defaults to 0: {json}");
        assert!(json.contains("\"rt\":0"), "rt defaults to 0: {json}");
        assert!(
            json.contains("\"n_points\":0"),
            "n_points (snake_case): {json}"
        );
        assert!(json.contains("\"noise\":0"), "noise field present: {json}");
        assert!(!json.contains("null"), "no null values: {json}");
        assert!(!json.contains("\"r2\""), "r2 is skipped: {json}");
    }

    #[test]
    fn test_spectrum_summary_unknown_serialises_without_null() {
        let s = SpectrumSummary::unknown();
        let json = serde_json::to_string(&s).expect("serialise ok");
        assert!(
            !json.contains("null"),
            "unknown SpectrumSummary must not produce null: {json}"
        );
        assert!(json.contains("\"rt_seconds\":0"), "rt_seconds → 0: {json}");
        assert!(
            json.contains("\"base_peak_mz\":0"),
            "base_peak_mz → 0: {json}"
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
            .map(|_| to_ion_bytes(&build_mzml_masses(&[mz_a, mz_b], peak_rt, 10_000.0, 4.0, 6.0, 60)))
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
            near, 1,
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
            assert_eq!(count_near(&features, mz, peak_rt), 1, "mass {mz} once, got {reported:?}");
        }
    }

    #[test]
    fn test_wider_tolerance_keeps_masses_merged() {
        let mz_a = 200.0300;
        let mz_b = 200.0330;
        let peak_rt = 5.0;
        let samples: Vec<Vec<u8>> = (0..6)
            .map(|_| to_ion_bytes(&build_mzml_masses(&[mz_a, mz_b], peak_rt, 10_000.0, 4.0, 6.0, 60)))
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
            near, 1,
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
}
