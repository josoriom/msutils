pub type Point = Vec<f64>;

/// K-means clustering (Lloyd’s algorithm).
/// Pseudocode reference: Wikipedia — “K-means clustering”
/// <https://en.wikipedia.org/wiki/K-means_clustering>

pub fn kmeans<P>(points: &[P], mut centroids: Vec<Point>) -> Vec<Point>
where
    P: AsRef<[f64]>,
{
    if points.is_empty() || centroids.is_empty() {
        return Vec::new();
    }

    let k = centroids.len();
    let d = centroids[0].len();
    if d == 0 {
        return Vec::new();
    }

    let mut converged = false;
    let mut iter = 0usize;

    let mut clusters: Vec<Vec<usize>> = (0..k).map(|_| Vec::new()).collect();
    let approx = points.len() / k + 1;
    for c in clusters.iter_mut() {
        c.reserve(approx);
    }

    let mut new_centroids: Vec<Point> = (0..k).map(|_| vec![0.0; d]).collect();

    let tol2 = 1e-12_f64;

    while !converged {
        for c in clusters.iter_mut() {
            c.clear();
        }

        for (i, p) in points.iter().enumerate() {
            let point = p.as_ref();

            let mut closest_index = 0usize;
            let mut min_distance = distance2(point, centroids[0].as_slice());

            for j in 1..k {
                let d2 = distance2(point, centroids[j].as_slice());
                if d2 < min_distance {
                    min_distance = d2;
                    closest_index = j;
                }
            }

            clusters[closest_index].push(i);
        }

        for i in 0..k {
            new_centroids[i].as_mut_slice().fill(0.0);

            if clusters[i].is_empty() {
                new_centroids[i]
                    .as_mut_slice()
                    .copy_from_slice(centroids[i].as_slice());
            } else {
                for &ix in clusters[i].iter() {
                    let p = points[ix].as_ref();
                    for t in 0..d {
                        new_centroids[i][t] += p[t];
                    }
                }

                let inv = 1.0 / clusters[i].len() as f64;
                for t in 0..d {
                    new_centroids[i][t] *= inv;
                }
            }
        }

        let mut max_shift2 = 0.0f64;
        for i in 0..k {
            let s2 = distance2(centroids[i].as_slice(), new_centroids[i].as_slice());
            if s2 > max_shift2 {
                max_shift2 = s2;
            }
        }

        converged = max_shift2 <= tol2;

        core::mem::swap(&mut centroids, &mut new_centroids);

        iter += 1;
        if iter > 300 {
            break;
        }
    }

    centroids
}

#[inline]
fn distance2(a: &[f64], b: &[f64]) -> f64 {
    let mut sum = 0.0;
    for i in 0..a.len() {
        let d = a[i] - b[i];
        sum += d * d;
    }
    sum
}

pub fn calculate_centroid(points_in_cluster: &[&Point]) -> Point {
    let d = points_in_cluster[0].len();
    let mut centroid = vec![0.0; d];

    for p in points_in_cluster {
        for i in 0..d {
            centroid[i] += p[i];
        }
    }
    let n = points_in_cluster.len() as f64;
    for i in 0..d {
        centroid[i] /= n;
    }
    centroid
}

pub fn distance(a: &Point, b: &Point) -> f64 {
    let mut sum = 0.0;
    for i in 0..a.len() {
        let d = a[i] - b[i];
        sum += d * d;
    }
    sum.sqrt()
}
