#[cfg(test)]
mod tests {
    use crate::utilities::{
        find_features::MzTolerance,
        mz_estimator::{
            MedianMzApex, MzEstimator, MzEstimatorKind, SampleMz, WeightedIntensityMedian,
            has_multiple_masses, highest_intensity, keep_near_anchor, make_estimator,
            same_mass_gap, split_by_gap, weighted_median,
        },
    };

    fn f06961_centroids() -> Vec<f64> {
        vec![
            119.02877, 119.02881, 119.02881, 119.02875, 119.02879, 119.02873, 119.02870,
            119.02885, 119.02894, 119.02895, 119.03149, 119.03152, 119.03153, 119.03150, 119.03153,
            119.03152, 119.03148, 119.03160, 119.03167, 119.03166, 119.02556, 119.02554,
        ]
    }

    fn f06961_gap() -> f64 {
        same_mass_gap(119.0288, 20.0)
    }

    fn group_median(group: &[f64]) -> f64 {
        let mut sorted = group.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        sorted[sorted.len() / 2]
    }

    #[test]
    fn test_has_multiple_masses_detects_f06961_merge() {
        assert!(has_multiple_masses(&f06961_centroids(), f06961_gap()));
    }

    #[test]
    fn test_split_by_gap_splits_f06961_two_masses() {
        let groups = split_by_gap(&f06961_centroids(), f06961_gap());
        let masses: Vec<&Vec<f64>> = groups.iter().filter(|g| g.len() >= 5).collect();
        assert_eq!(masses.len(), 2, "exactly two real masses");

        let mut medians: Vec<f64> = masses.iter().map(|g| group_median(g)).collect();
        medians.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((medians[0] - 119.0288).abs() < 0.0005, "low mass {}", medians[0]);
        assert!((medians[1] - 119.0315).abs() < 0.0005, "high mass {}", medians[1]);

        let merged = groups
            .iter()
            .any(|g| g.iter().any(|&m| m < 119.030) && g.iter().any(|&m| m > 119.031));
        assert!(!merged, "the two masses must never share one group");
    }

    #[test]
    fn test_split_by_gap_isolates_outlier() {
        let groups = split_by_gap(&f06961_centroids(), f06961_gap());
        assert!(
            groups
                .iter()
                .any(|g| g.len() == 2 && group_median(g) < 119.027),
            "the 119.0255 outlier is its own small group"
        );
    }

    #[test]
    fn test_has_multiple_masses_empty_is_false() {
        assert!(!has_multiple_masses(&[], f06961_gap()));
        assert!(split_by_gap(&[], f06961_gap()).is_empty());
    }

    #[test]
    fn test_has_multiple_masses_single_is_false() {
        assert!(!has_multiple_masses(&[119.0288], f06961_gap()));
        assert_eq!(split_by_gap(&[119.0288], f06961_gap()).len(), 1);
    }

    #[test]
    fn test_one_tight_cluster_not_split() {
        let centroids = [119.02870, 119.02875, 119.02880, 119.02885, 119.02890];
        assert!(!has_multiple_masses(&centroids, f06961_gap()));
        assert_eq!(split_by_gap(&centroids, f06961_gap()).len(), 1);
    }

    #[test]
    fn test_two_masses_within_cutoff_stay_merged() {
        let centroids = [119.0288, 119.0289, 119.0296, 119.0297];
        assert!(!has_multiple_masses(&centroids, f06961_gap()));
        assert_eq!(split_by_gap(&centroids, f06961_gap()).len(), 1);
    }

    #[test]
    fn test_two_clear_masses_split() {
        let centroids = [119.0287, 119.0288, 119.0289, 119.0314, 119.0315, 119.0316];
        assert!(has_multiple_masses(&centroids, f06961_gap()));
        assert_eq!(split_by_gap(&centroids, f06961_gap()).len(), 2);
    }

    #[test]
    fn test_three_masses_split() {
        let centroids = [100.0000, 100.0001, 100.0050, 100.0051, 100.0100, 100.0101];
        let gap = same_mass_gap(100.0, 20.0);
        assert_eq!(split_by_gap(&centroids, gap).len(), 3);
    }

    #[test]
    fn test_same_mass_gap_is_ppm_based() {
        assert!((same_mass_gap(100.0, 20.0) - 0.001).abs() < 1e-12);
        assert!((same_mass_gap(200.0, 20.0) - 0.002).abs() < 1e-12);
        assert!(same_mass_gap(119.0, 10.0) < same_mass_gap(119.0, 20.0));
    }

    fn tol() -> MzTolerance {
        MzTolerance {
            mz_absolute: 0.005,
            ppm: 0.0,
        }
    }

    fn sample(mz: f64, intensity: f64) -> SampleMz {
        SampleMz { mz, intensity }
    }

    fn scan(rt: f64, mz: Vec<f64>, intensity: Vec<f64>) -> (f64, Vec<f64>, Vec<f64>) {
        (rt, mz, intensity)
    }

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-9
    }

    #[test]
    fn test_weighted_median_empty_is_none() {
        assert!(weighted_median(&[]).is_none());
    }

    #[test]
    fn test_weighted_median_single() {
        assert!(close(weighted_median(&[sample(42.0, 3.0)]).unwrap(), 42.0));
    }

    #[test]
    fn test_weighted_median_three_equal_weights_picks_middle() {
        let values = [sample(10.0, 1.0), sample(20.0, 1.0), sample(30.0, 1.0)];
        assert!(close(weighted_median(&values).unwrap(), 20.0));
    }

    #[test]
    fn test_weighted_median_crossing_exactly_at_half() {
        let values = [sample(10.0, 1.0), sample(20.0, 1.0)];
        assert!(close(weighted_median(&values).unwrap(), 10.0));
    }

    #[test]
    fn test_weighted_median_dominant_weight_wins() {
        let values = [sample(10.0, 0.0), sample(20.0, 5.0), sample(30.0, 0.0)];
        assert!(close(weighted_median(&values).unwrap(), 20.0));
    }

    #[test]
    fn test_weighted_median_zero_total_is_none() {
        let values = [sample(10.0, 0.0), sample(20.0, 0.0)];
        assert!(weighted_median(&values).is_none());
    }

    #[test]
    fn test_weighted_median_unsorted_input() {
        let values = [sample(30.0, 1.0), sample(10.0, 1.0), sample(20.0, 1.0)];
        assert!(close(weighted_median(&values).unwrap(), 20.0));
    }

    #[test]
    fn test_highest_intensity_empty_is_none() {
        assert!(highest_intensity(&[]).is_none());
    }

    #[test]
    fn test_highest_intensity_picks_max() {
        let values = [sample(10.0, 5.0), sample(20.0, 9.0), sample(30.0, 1.0)];
        let top = highest_intensity(&values).unwrap();
        assert!(close(top.mz, 20.0));
        assert!(close(top.intensity, 9.0));
    }

    #[test]
    fn test_highest_intensity_tie_keeps_one_with_max() {
        let values = [sample(10.0, 5.0), sample(20.0, 5.0)];
        assert!(close(highest_intensity(&values).unwrap().intensity, 5.0));
    }

    #[test]
    fn test_keep_near_anchor_empty_is_empty() {
        assert!(keep_near_anchor(&[], &tol()).is_empty());
    }

    #[test]
    fn test_keep_near_anchor_drops_outliers() {
        let values = [
            sample(100.000, 50.0),
            sample(100.002, 90.0),
            sample(100.004, 40.0),
            sample(100.500, 10.0),
        ];
        let kept = keep_near_anchor(&values, &tol());
        assert_eq!(kept.len(), 3);
        assert!(kept.iter().all(|v| v.mz < 100.1));
    }

    #[test]
    fn test_keep_near_anchor_all_within_tolerance() {
        let values = [sample(100.000, 1.0), sample(100.001, 2.0)];
        assert_eq!(keep_near_anchor(&values, &tol()).len(), 2);
    }

    #[test]
    fn test_apex_per_sample_picks_apex_scan_most_intense_point() {
        let scans = [
            scan(4.9, vec![100.000, 100.004], vec![10.0, 5.0]),
            scan(5.0, vec![100.001, 100.005], vec![80.0, 20.0]),
            scan(5.1, vec![100.002, 100.006], vec![15.0, 7.0]),
        ];
        let result = MedianMzApex.per_sample(&scans, 99.99, 100.01).unwrap();
        assert!(close(result.mz, 100.001));
        assert!(close(result.intensity, 80.0));
    }

    #[test]
    fn test_apex_per_sample_single_scan() {
        let scans = [scan(5.0, vec![100.001, 100.005], vec![80.0, 20.0])];
        let result = MedianMzApex.per_sample(&scans, 99.99, 100.01).unwrap();
        assert!(close(result.mz, 100.001));
    }

    #[test]
    fn test_apex_per_sample_includes_point_on_window_edge() {
        let scans = [scan(5.0, vec![100.001, 100.005], vec![20.0, 80.0])];
        let result = MedianMzApex.per_sample(&scans, 99.99, 100.005).unwrap();
        assert!(close(result.mz, 100.005));
    }

    #[test]
    fn test_apex_per_sample_no_points_in_window_is_none() {
        let scans = [scan(5.0, vec![200.0, 201.0], vec![80.0, 20.0])];
        assert!(MedianMzApex.per_sample(&scans, 99.99, 100.01).is_none());
    }

    #[test]
    fn test_apex_per_sample_all_zero_intensity_is_none() {
        let scans = [scan(5.0, vec![100.001, 100.002], vec![0.0, 0.0])];
        assert!(MedianMzApex.per_sample(&scans, 99.99, 100.01).is_none());
    }

    #[test]
    fn test_apex_per_sample_skips_non_finite_points() {
        let scans = [scan(
            5.0,
            vec![f64::NAN, 100.002],
            vec![999.0, 50.0],
        )];
        let result = MedianMzApex.per_sample(&scans, 99.99, 100.01).unwrap();
        assert!(close(result.mz, 100.002));
        assert!(close(result.intensity, 50.0));
    }

    #[test]
    fn test_weighted_median_non_finite_total_is_none() {
        let values = [sample(100.0, f64::INFINITY)];
        assert!(weighted_median(&values).is_none());
    }

    #[test]
    fn test_apex_combine_drops_outlier_and_takes_median() {
        let values = [
            sample(100.000, 50.0),
            sample(100.002, 90.0),
            sample(100.004, 40.0),
            sample(100.500, 10.0),
        ];
        let result = MedianMzApex.combine(&values, &tol()).unwrap();
        assert!(close(result, 100.002));
    }

    #[test]
    fn test_apex_combine_even_count_averages_middle() {
        let values = [sample(100.000, 30.0), sample(100.002, 90.0)];
        let result = MedianMzApex.combine(&values, &tol()).unwrap();
        assert!(close(result, 100.001));
    }

    #[test]
    fn test_apex_combine_single() {
        let result = MedianMzApex.combine(&[sample(100.002, 90.0)], &tol()).unwrap();
        assert!(close(result, 100.002));
    }

    #[test]
    fn test_apex_combine_empty_is_none() {
        assert!(MedianMzApex.combine(&[], &tol()).is_none());
    }

    #[test]
    fn test_apex_combine_anchor_follows_most_intense_even_when_isolated() {
        let values = [
            sample(100.000, 10.0),
            sample(100.002, 12.0),
            sample(100.004, 11.0),
            sample(100.030, 500.0),
        ];
        let result = MedianMzApex.combine(&values, &tol()).unwrap();
        assert!(close(result, 100.030));
    }

    #[test]
    fn test_weighted_per_sample_intensity_weighted_median_over_window() {
        let scans = [
            scan(5.0, vec![100.000, 100.004], vec![1.0, 1.0]),
            scan(5.1, vec![100.008], vec![10.0]),
        ];
        let result = WeightedIntensityMedian
            .per_sample(&scans, 99.99, 100.01)
            .unwrap();
        assert!(close(result.mz, 100.008));
        assert!(close(result.intensity, 12.0));
    }

    #[test]
    fn test_weighted_per_sample_single_point() {
        let scans = [scan(5.0, vec![100.002], vec![5.0])];
        let result = WeightedIntensityMedian
            .per_sample(&scans, 99.99, 100.01)
            .unwrap();
        assert!(close(result.mz, 100.002));
    }

    #[test]
    fn test_weighted_per_sample_zero_intensity_is_none() {
        let scans = [scan(5.0, vec![100.002], vec![0.0])];
        assert!(
            WeightedIntensityMedian
                .per_sample(&scans, 99.99, 100.01)
                .is_none()
        );
    }

    #[test]
    fn test_weighted_combine_drops_outlier() {
        let values = [
            sample(100.000, 1.0),
            sample(100.002, 9.0),
            sample(100.500, 2.0),
        ];
        let result = WeightedIntensityMedian.combine(&values, &tol()).unwrap();
        assert!(close(result, 100.002));
    }

    #[test]
    fn test_weighted_combine_empty_is_none() {
        assert!(WeightedIntensityMedian.combine(&[], &tol()).is_none());
    }

    #[test]
    fn test_make_estimator_median_apex_matches_direct() {
        let estimator = make_estimator(&MzEstimatorKind::MedianMzApex);
        let values = [sample(100.000, 50.0), sample(100.002, 90.0)];
        let from_kind = estimator.combine(&values, &tol()).unwrap();
        let direct = MedianMzApex.combine(&values, &tol()).unwrap();
        assert!(close(from_kind, direct));
    }

    #[test]
    fn test_make_estimator_weighted_matches_direct() {
        let estimator = make_estimator(&MzEstimatorKind::WeightedIntensityMedian);
        let values = [sample(100.000, 1.0), sample(100.002, 9.0)];
        let from_kind = estimator.combine(&values, &tol()).unwrap();
        let direct = WeightedIntensityMedian.combine(&values, &tol()).unwrap();
        assert!(close(from_kind, direct));
    }
}
