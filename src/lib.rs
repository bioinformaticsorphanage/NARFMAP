//! NARFMAP - DRAGEN-OS mapper/aligner library
//!
//! This library provides genomic sequence alignment functionality.

pub mod alignment;
pub mod ffi;
pub mod io;
pub mod reference;

/// Alignment scoring parameters
#[derive(Debug, Clone)]
pub struct ScoringParams {
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

impl Default for ScoringParams {
    fn default() -> Self {
        Self {
            match_score: 1,
            mismatch_score: -4,
            gap_open: -7,   // bio crate expects non-positive penalty
            gap_extend: -1, // bio crate expects non-positive penalty
        }
    }
}

/// A genomic read with sequence and quality scores
#[derive(Debug, Clone)]
pub struct Read {
    pub name: String,
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
}

/// An alignment result
#[derive(Debug, Clone)]
pub struct Alignment {
    pub read_name: String,
    pub flag: u16,
    pub ref_name: String,
    pub position: u64,
    pub mapq: u8,
    pub cigar: String,
    pub sequence: Vec<u8>,
    pub quality: Vec<u8>,
    pub score: i32,
}

impl Alignment {
    /// Create an unmapped alignment
    pub fn unmapped(read: &Read) -> Self {
        Self {
            read_name: read.name.clone(),
            flag: 4, // unmapped
            ref_name: "*".to_string(),
            position: 0,
            mapq: 0,
            cigar: "*".to_string(),
            sequence: read.sequence.clone(),
            quality: read.quality.clone(),
            score: 0,
        }
    }
}
