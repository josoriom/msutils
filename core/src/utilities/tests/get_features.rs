#[cfg(test)]
mod tests {
    use crate::utilities::{
        find_features::{Feature, MzTolerance},
        get_features::{
            FeatureClusterer, GrowingCluster, SearchBounds, TaggedFeature,
            aggregate_into_consensus, assign_best_per_sample, collect_filled_slots,
            compute_search_bounds, median, require_minimum_prevalence, rsd,
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

    fn make_feature_full(mz: f64, rt: f64, intensity: f64, np: usize, integral: f64) -> Feature {
        Feature {
            mz,
            rt,
            intensity,
            np,
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
                mz_abs: 0.01,
                ppm: 0.0,
            },
            rt_tolerance: 0.1,
        }
    }

    fn default_bounds() -> SearchBounds {
        SearchBounds {
            target_mz: 100.0,
            seed_rt: 1.0,
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
    fn test_rsd_uniform_values() {
        assert_eq!(rsd(&[5.0, 5.0, 5.0]), 0.0);
    }

    #[test]
    fn test_rsd_single_value() {
        assert_eq!(rsd(&[5.0]), 0.0);
    }

    #[test]
    fn test_rsd_zero_mean() {
        assert_eq!(rsd(&[0.0, 0.0]), 0.0);
    }

    #[test]
    fn test_rsd_known_values() {
        let v = rsd(&[1.0, 2.0, 3.0]);
        assert!((v - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_rsd_two_elements() {
        let v = rsd(&[1.0, 2.0]);
        let expected = 0.5f64.sqrt() / 1.5;
        assert!((v - expected).abs() < 1e-10);
    }

    #[test]
    fn test_growing_cluster_single_item() {
        let gc = GrowingCluster::new(make_tagged(0, 100.0, 1.0, 500.0));
        assert_eq!(gc.cached_median_mz, 100.0);
        assert_eq!(gc.cached_median_rt, 1.0);
    }

    #[test]
    fn test_growing_cluster_median_updates_on_push() {
        let mut gc = GrowingCluster::new(make_tagged(0, 100.0, 1.0, 500.0));
        gc.push(make_tagged(1, 102.0, 2.0, 300.0));
        assert_eq!(gc.cached_median_mz, 102.0);
    }

    #[test]
    fn test_growing_cluster_three_items_median() {
        let mut gc = GrowingCluster::new(make_tagged(0, 100.0, 1.0, 500.0));
        gc.push(make_tagged(1, 102.0, 2.0, 300.0));
        gc.push(make_tagged(2, 104.0, 3.0, 200.0));
        assert_eq!(gc.cached_median_mz, 102.0);
    }

    #[test]
    fn test_growing_cluster_four_items_upper_index() {
        let mut gc = GrowingCluster::new(make_tagged(0, 100.0, 1.0, 500.0));
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
        assert_eq!(bounds.seed_rt, 1.0);
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
        assert_eq!(bounds.seed_rt, 2.0);
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
    fn test_require_minimum_prevalence_exact_threshold() {
        let hits = vec![make_feature(100.0, 1.0, 500.0)];
        assert!(require_minimum_prevalence(hits, 1).is_some());
    }

    #[test]
    fn test_require_minimum_prevalence_passes() {
        let hits = vec![
            make_feature(100.0, 1.0, 500.0),
            make_feature(100.0, 1.0, 400.0),
        ];
        assert!(require_minimum_prevalence(hits, 2).is_some());
    }

    #[test]
    fn test_require_minimum_prevalence_fails() {
        let hits = vec![make_feature(100.0, 1.0, 500.0)];
        assert!(require_minimum_prevalence(hits, 3).is_none());
    }

    #[test]
    fn test_require_minimum_prevalence_empty() {
        assert!(require_minimum_prevalence(vec![], 1).is_none());
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
        let result = aggregate_into_consensus(hits, &default_bounds());
        assert_eq!(result.frequency, 3);
    }

    #[test]
    fn test_aggregate_from_to_from_bounds() {
        let hits = vec![make_feature(100.0, 1.0, 500.0)];
        let bounds = SearchBounds {
            target_mz: 100.0,
            seed_rt: 1.0,
            rt_from: 0.8,
            rt_to: 1.2,
        };
        let result = aggregate_into_consensus(hits, &bounds);
        assert_eq!(result.from, 0.8);
        assert_eq!(result.to, 1.2);
    }

    #[test]
    fn test_aggregate_intensity_even_median() {
        let hits = vec![
            make_feature(100.0, 1.0, 300.0),
            make_feature(100.0, 1.0, 500.0),
        ];
        let result = aggregate_into_consensus(hits, &default_bounds());
        assert_eq!(result.intensity, 400.0);
    }

    #[test]
    fn test_aggregate_mz_odd_median() {
        let hits = vec![
            make_feature(99.0, 1.0, 500.0),
            make_feature(100.0, 1.0, 400.0),
            make_feature(101.0, 1.0, 300.0),
        ];
        let result = aggregate_into_consensus(hits, &default_bounds());
        assert_eq!(result.mz, 100.0);
    }

    #[test]
    fn test_aggregate_np_median() {
        let hits = vec![
            make_feature_full(100.0, 1.0, 500.0, 2, 10.0),
            make_feature_full(100.0, 1.0, 400.0, 4, 20.0),
            make_feature_full(100.0, 1.0, 300.0, 6, 30.0),
        ];
        let result = aggregate_into_consensus(hits, &default_bounds());
        assert_eq!(result.np, 4);
    }

    #[test]
    fn test_aggregate_integral_median() {
        let hits = vec![
            make_feature_full(100.0, 1.0, 500.0, 2, 10.0),
            make_feature_full(100.0, 1.0, 400.0, 4, 20.0),
            make_feature_full(100.0, 1.0, 300.0, 6, 30.0),
        ];
        let result = aggregate_into_consensus(hits, &default_bounds());
        assert!((result.integral - 20.0).abs() < 1e-10);
    }

    #[test]
    fn test_aggregate_correct_frequency() {
        let hits = vec![
            make_feature(100.0, 1.0, 500.0),
            make_feature(100.01, 1.05, 400.0),
            make_feature(99.99, 0.95, 600.0),
        ];
        let result = aggregate_into_consensus(hits, &default_bounds());
        assert_eq!(result.frequency, 3);
        assert!((result.mz - 100.0).abs() < 0.01);
    }
}
