use ionic::ion::ByteRange;

use crate::utilities::calculate_eic::FastError;

pub type RangePlan<'a> = &'a mut dyn FnMut() -> Result<Vec<ByteRange>, FastError>;

pub trait Prefetcher: Send + Sync {
    fn prefetch(&self, plan: RangePlan) -> Result<(), FastError>;
}

pub struct NoPrefetch;

impl Prefetcher for NoPrefetch {
    fn prefetch(&self, _plan: RangePlan) -> Result<(), FastError> {
        Ok(())
    }
}
