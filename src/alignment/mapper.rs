//! Seed-based read mapping

use std::collections::HashMap;

use crate::reference::hash_table::{HashTable, RefPosition};
use crate::Read;

/// A candidate alignment region
#[derive(Debug, Clone)]
pub struct CandidateRegion {
    pub ref_start: u64,
    pub ref_end: u64,
    pub seed_count: usize,
    pub is_reverse: bool,
}

/// Seed-based mapper for finding candidate alignment regions
pub struct SeedMapper<'a> {
    hash_table: &'a HashTable,
    seed_step: usize,
    min_seeds: usize,
}

impl<'a> SeedMapper<'a> {
    /// Create a new seed mapper
    pub fn new(hash_table: &'a HashTable) -> Self {
        Self {
            hash_table,
            seed_step: 1, // Step between seeds
            min_seeds: 2, // Minimum seeds to form a candidate
        }
    }

    /// Set the step size between seeds
    pub fn with_seed_step(mut self, step: usize) -> Self {
        self.seed_step = step.max(1);
        self
    }

    /// Set minimum seeds required for a candidate region
    pub fn with_min_seeds(mut self, min: usize) -> Self {
        self.min_seeds = min.max(1);
        self
    }

    /// Map a read to candidate regions
    pub fn map(&self, read: &Read) -> Vec<CandidateRegion> {
        let seed_len = self.hash_table.config.seed_len();

        if read.sequence.len() < seed_len {
            return Vec::new();
        }

        // Collect all seed hits
        let mut hits_fwd: HashMap<u64, Vec<(usize, u64)>> = HashMap::new(); // bin -> (read_pos, ref_pos)
        let mut hits_rev: HashMap<u64, Vec<(usize, u64)>> = HashMap::new();

        let bin_size = 100u64; // Group hits within 100bp bins

        for read_pos in (0..=read.sequence.len() - seed_len).step_by(self.seed_step) {
            let seed = &read.sequence[read_pos..read_pos + seed_len];

            for hit in self.hash_table.query(seed) {
                let bin = hit.offset / bin_size;

                if hit.is_reverse {
                    hits_rev
                        .entry(bin)
                        .or_default()
                        .push((read_pos, hit.offset));
                } else {
                    hits_fwd
                        .entry(bin)
                        .or_default()
                        .push((read_pos, hit.offset));
                }
            }
        }

        // Find candidate regions from clustered hits
        let mut candidates = Vec::new();

        for (_, hits) in hits_fwd {
            if hits.len() >= self.min_seeds {
                let ref_start = hits.iter().map(|(_, p)| *p).min().unwrap_or(0);
                let ref_end = hits.iter().map(|(_, p)| *p).max().unwrap_or(0) + seed_len as u64;

                candidates.push(CandidateRegion {
                    ref_start,
                    ref_end,
                    seed_count: hits.len(),
                    is_reverse: false,
                });
            }
        }

        for (_, hits) in hits_rev {
            if hits.len() >= self.min_seeds {
                let ref_start = hits.iter().map(|(_, p)| *p).min().unwrap_or(0);
                let ref_end = hits.iter().map(|(_, p)| *p).max().unwrap_or(0) + seed_len as u64;

                candidates.push(CandidateRegion {
                    ref_start,
                    ref_end,
                    seed_count: hits.len(),
                    is_reverse: true,
                });
            }
        }

        // Sort by seed count (descending)
        candidates.sort_by(|a, b| b.seed_count.cmp(&a.seed_count));

        // Limit to top candidates
        candidates.truncate(10);

        candidates
    }
}

/// Reverse complement a DNA sequence
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            _ => b'N',
        })
        .collect()
}
