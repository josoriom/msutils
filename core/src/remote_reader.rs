use std::sync::Mutex;

use ionic::ion::{ByteRange, IonResult, open_ranges};
use rayon::prelude::*;

use crate::prefetch::{Prefetcher, RangePlan};
use crate::url_source::UrlSource;
use crate::utilities::calculate_eic::{FastError, sort_and_dedup_ranges};

const HEADER_LENGTH: u64 = 1024;
const MERGE_GAP: u64 = 65536;

struct CachedBytes {
    start: u64,
    end: u64,
    bytes: Vec<u8>,
}

struct RangeCache {
    stored: Mutex<Vec<CachedBytes>>,
}

impl RangeCache {
    fn new() -> Self {
        Self {
            stored: Mutex::new(Vec::new()),
        }
    }

    fn get(&self, offset: u64, length: u64) -> Option<Vec<u8>> {
        let end = offset + length;
        let stored = self.stored.lock().unwrap();
        let count = stored.partition_point(|range| range.start <= offset);
        if count == 0 {
            return None;
        }
        let range = &stored[count - 1];
        if range.end < end {
            return None;
        }
        let from = (offset - range.start) as usize;
        let to = (end - range.start) as usize;
        Some(range.bytes[from..to].to_vec())
    }

    fn add(&self, offset: u64, bytes: Vec<u8>) {
        let entry = CachedBytes {
            start: offset,
            end: offset + bytes.len() as u64,
            bytes,
        };
        let mut stored = self.stored.lock().unwrap();
        let at = stored.partition_point(|range| range.start < entry.start);
        debug_assert!(
            at == 0 || stored[at - 1].end <= entry.start,
            "cached ranges must stay non-overlapping"
        );
        debug_assert!(
            at == stored.len() || entry.end <= stored[at].start,
            "cached ranges must stay non-overlapping"
        );
        stored.insert(at, entry);
    }

    fn missing(&self, wanted: &[ByteRange]) -> Vec<ByteRange> {
        wanted
            .iter()
            .filter(|range| self.get(range.offset, range.length).is_none())
            .copied()
            .collect()
    }

    fn clear(&self) {
        self.stored.lock().unwrap().clear();
    }
}

fn merge_ranges(ranges: &[ByteRange]) -> Vec<ByteRange> {
    if ranges.is_empty() {
        return Vec::new();
    }
    let mut sorted = ranges.to_vec();
    sort_and_dedup_ranges(&mut sorted);
    let mut merged = Vec::with_capacity(sorted.len());
    let mut current = sorted[0];
    for next in sorted.into_iter().skip(1) {
        let current_end = current.offset + current.length;
        if next.offset <= current_end + MERGE_GAP {
            let next_end = next.offset + next.length;
            current.length = current_end.max(next_end) - current.offset;
        } else {
            merged.push(current);
            current = next;
        }
    }
    merged.push(current);
    merged
}

pub struct RemoteReader {
    source: UrlSource,
    cache: RangeCache,
}

impl RemoteReader {
    pub fn new(url: String) -> IonResult<Self> {
        Ok(Self {
            source: UrlSource::new(url)?,
            cache: RangeCache::new(),
        })
    }

    pub fn read(&self, range: ByteRange) -> IonResult<Vec<u8>> {
        let length = range.length;
        if length == 0 {
            return Ok(Vec::new());
        }
        if let Some(bytes) = self.cache.get(range.offset, length) {
            return Ok(bytes);
        }
        self.source.read(range)
    }

    pub fn prefetch_open(&self) -> IonResult<()> {
        let header = self.source.read(ByteRange { offset: 0, length: HEADER_LENGTH })?;
        let ranges = open_ranges(&header)?;
        self.cache.add(0, header);
        self.fetch_into_cache(&ranges)
    }

    pub fn clear(&self) {
        self.cache.clear();
    }

    fn fetch_into_cache(&self, ranges: &[ByteRange]) -> IonResult<()> {
        let merged = merge_ranges(&self.cache.missing(ranges));
        let fetched: IonResult<Vec<(u64, Vec<u8>)>> = merged
            .par_iter()
            .map(|range| {
                let bytes = self.source.read(ByteRange { offset: range.offset, length: range.length })?;
                Ok((range.offset, bytes))
            })
            .collect();
        for (offset, bytes) in fetched? {
            self.cache.add(offset, bytes);
        }
        Ok(())
    }
}

impl Prefetcher for RemoteReader {
    fn prefetch(&self, plan: RangePlan) -> Result<(), FastError> {
        let ranges = plan()?;
        self.cache.clear();
        self.fetch_into_cache(&ranges)
            .map_err(|error| FastError::ReadFailed(error.to_string()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn range(offset: u64, length: u64) -> ByteRange {
        ByteRange { offset, length }
    }

    #[test]
    fn cache_returns_exact_bytes() {
        let cache = RangeCache::new();
        cache.add(100, vec![1, 2, 3, 4]);
        assert_eq!(cache.get(100, 4), Some(vec![1, 2, 3, 4]));
    }

    #[test]
    fn cache_returns_slice_inside_larger_range() {
        let cache = RangeCache::new();
        cache.add(100, vec![1, 2, 3, 4, 5, 6]);
        assert_eq!(cache.get(102, 3), Some(vec![3, 4, 5]));
    }

    #[test]
    fn cache_misses_when_nothing_covers_request() {
        let cache = RangeCache::new();
        cache.add(100, vec![1, 2, 3, 4]);
        assert_eq!(cache.get(200, 4), None);
    }

    #[test]
    fn cache_misses_when_request_runs_past_end() {
        let cache = RangeCache::new();
        cache.add(100, vec![1, 2, 3, 4]);
        assert_eq!(cache.get(102, 4), None);
    }

    #[test]
    fn cache_keeps_ranges_sorted_by_start() {
        let cache = RangeCache::new();
        cache.add(300, vec![3]);
        cache.add(100, vec![1]);
        cache.add(200, vec![2]);
        assert_eq!(cache.get(100, 1), Some(vec![1]));
        assert_eq!(cache.get(200, 1), Some(vec![2]));
        assert_eq!(cache.get(300, 1), Some(vec![3]));
    }

    #[test]
    fn missing_keeps_only_uncovered_ranges() {
        let cache = RangeCache::new();
        cache.add(100, vec![0; 50]);
        let wanted = vec![range(100, 10), range(500, 10)];
        assert_eq!(cache.missing(&wanted), vec![range(500, 10)]);
    }

    #[test]
    fn merge_joins_ranges_within_gap() {
        let merged = merge_ranges(&[range(0, 100), range(100 + MERGE_GAP, 100)]);
        assert_eq!(merged, vec![range(0, 100 + MERGE_GAP + 100)]);
    }

    #[test]
    fn merge_keeps_distant_ranges_apart() {
        let far = 100 + MERGE_GAP + 1;
        let merged = merge_ranges(&[range(0, 100), range(far, 100)]);
        assert_eq!(merged, vec![range(0, 100), range(far, 100)]);
    }

    #[test]
    fn merge_collapses_overlapping_ranges() {
        let merged = merge_ranges(&[range(0, 200), range(150, 100)]);
        assert_eq!(merged, vec![range(0, 250)]);
    }

    #[test]
    fn merge_handles_empty_input() {
        assert_eq!(merge_ranges(&[]), Vec::new());
    }
}
