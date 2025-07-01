use std::fmt;
use anyhow::Result;
// Use our stub types from reference module instead of bitnuc
use crate::reference::sequence::{NucSeq, Nucleotide};

/// Represents a DNA sequence with an identifier and optional quality scores
pub struct Sequence {
    /// The sequence identifier (e.g., read name)
    pub id: String,
    /// The nucleotide sequence represented efficiently using bitnuc's 2-bit encoding
    pub seq: NucSeq,
    /// Optional quality scores
    pub qual: Option<Vec<u8>>,
}

impl Sequence {
    /// Create a new sequence from a string representation
    pub fn new(id: String, seq_str: &str, qual: Option<Vec<u8>>) -> Result<Self> {
        let seq = NucSeq::from_str(seq_str)?;
        Ok(Self { id, seq, qual })
    }

    /// Get the sequence length
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    /// Check if the sequence is empty
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    /// Get the nucleotide at the specified position
    pub fn get(&self, idx: usize) -> Option<Nucleotide> {
        if idx < self.len() {
            Some(self.seq.get(idx))
        } else {
            None
        }
    }

    /// Extract a subsequence (substring)
    pub fn substring(&self, start: usize, end: usize) -> Option<NucSeq> {
        if start <= end && end <= self.len() {
            Some(self.seq.subseq(start, end))
        } else {
            None
        }
    }

    /// Reverse complement the sequence
    pub fn reverse_complement(&self) -> Self {
        Self {
            id: self.id.clone(),
            seq: self.seq.clone(),
            qual: self.qual.as_ref().map(|q| q.iter().rev().cloned().collect()),
        }
    }
}

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, ">{}\n{}", self.id, self.seq)
    }
}

impl fmt::Debug for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Sequence")
            .field("id", &self.id)
            .field("seq", &self.seq.to_string())
            .field("qual", &self.qual.as_ref().map(|q| format!("{:?}", q)))
            .finish()
    }
} 