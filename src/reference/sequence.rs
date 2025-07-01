use std::fmt;
use std::ops::Range;
use std::collections::HashMap;
use anyhow::{Result, Context, anyhow};
use bitnuc;

use crate::io::sequence::Sequence as IoSequence;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Nucleotide {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
    N = 4, // Unknown/ambiguous nucleotide
}

impl Nucleotide {
    pub fn from_byte(b: u8) -> Self {
        match b.to_ascii_uppercase() {
            b'A' => Nucleotide::A,
            b'C' => Nucleotide::C,
            b'G' => Nucleotide::G,
            b'T' => Nucleotide::T,
            _ => Nucleotide::N,
        }
    }

    pub fn to_byte(self) -> u8 {
        match self {
            Nucleotide::A => b'A',
            Nucleotide::C => b'C',
            Nucleotide::G => b'G',
            Nucleotide::T => b'T',
            Nucleotide::N => b'N',
        }
    }

    pub fn complement(self) -> Self {
        match self {
            Nucleotide::A => Nucleotide::T,
            Nucleotide::T => Nucleotide::A,
            Nucleotide::C => Nucleotide::G,
            Nucleotide::G => Nucleotide::C,
            Nucleotide::N => Nucleotide::N,
        }
    }
}

/// Efficient nucleotide sequence using bitnuc's 2-bit encoding
#[derive(Debug, Clone)]
pub struct NucSeq {
    /// Raw sequence data as bytes (for sequences with N's or other ambiguous bases)
    raw_seq: Vec<u8>,
    /// Length of the sequence
    length: usize,
}

impl NucSeq {
    pub fn from_str(s: &str) -> Result<Self> {
        let bytes = s.as_bytes().to_vec();
        Ok(Self {
            length: bytes.len(),
            raw_seq: bytes,
        })
    }

    pub fn from_bytes(bytes: &[u8]) -> Self {
        Self {
            length: bytes.len(),
            raw_seq: bytes.to_vec(),
        }
    }

    pub fn len(&self) -> usize {
        self.length
    }

    pub fn is_empty(&self) -> bool {
        self.length == 0
    }

    pub fn get(&self, pos: usize) -> Nucleotide {
        if pos < self.length {
            Nucleotide::from_byte(self.raw_seq[pos])
        } else {
            Nucleotide::N
        }
    }

    pub fn subseq(&self, start: usize, end: usize) -> Self {
        if start <= end && end <= self.length {
            Self {
                raw_seq: self.raw_seq[start..end].to_vec(),
                length: end - start,
            }
        } else {
            Self {
                raw_seq: Vec::new(),
                length: 0,
            }
        }
    }

    pub fn to_string(&self) -> String {
        String::from_utf8_lossy(&self.raw_seq).to_string()
    }

    /// Get raw bytes for the sequence
    pub fn as_bytes(&self) -> &[u8] {
        &self.raw_seq
    }

    /// Convert to 2-bit encoding for k-mers (up to 32 bases)
    pub fn to_2bit_kmer(&self, k: usize) -> Result<u64> {
        if k > 32 {
            return Err(anyhow!("K-mer size {} too large for 2-bit encoding (max 32)", k));
        }
        
        if self.length < k {
            return Err(anyhow!("Sequence length {} shorter than k-mer size {}", self.length, k));
        }

        // For sequences with N's, we need to handle them differently
        // For now, convert N's to A's for 2-bit encoding
        let clean_seq: Vec<u8> = self.raw_seq[0..k].iter()
            .map(|&b| match b.to_ascii_uppercase() {
                b'A' | b'C' | b'G' | b'T' => b,
                _ => b'A', // Convert N's and other ambiguous bases to A
            })
            .collect();

        bitnuc::as_2bit(&clean_seq).map_err(|e| anyhow!("Failed to encode as 2-bit: {}", e))
    }

    /// Reverse complement of the sequence
    pub fn reverse_complement(&self) -> Self {
        let rev_comp: Vec<u8> = self.raw_seq
            .iter()
            .rev()
            .map(|&b| Nucleotide::from_byte(b).complement().to_byte())
            .collect();
        
        Self {
            raw_seq: rev_comp,
            length: self.length,
        }
    }
}

impl fmt::Display for NucSeq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_string())
    }
}

/// Represents a reference sequence from a genome assembly
#[derive(Debug, Clone)]
pub struct ReferenceSequence {
    /// The sequence data stored efficiently using bitnuc's 2-bit encoding
    sequence: NucSeq,
    /// The name of the sequence (typically chromosome name)
    name: String,
    /// Length of the sequence
    length: usize,
}

impl ReferenceSequence {
    /// Create a new reference sequence from a bitnuc NucSeq
    pub fn new(name: &str, sequence: NucSeq) -> Self {
        let length = sequence.len();
        Self {
            sequence,
            name: name.to_string(),
            length,
        }
    }
    
    /// Create a new reference sequence from raw byte data
    pub fn from_bytes(name: &str, sequence: &[u8]) -> Result<Self> {
        // Convert raw bytes to a string then to a NucSeq
        let seq_str = std::str::from_utf8(sequence)
            .with_context(|| "Invalid UTF-8 sequence data")?;
        
        let nuc_seq = NucSeq::from_str(seq_str)
            .with_context(|| "Failed to convert sequence to NucSeq")?;
        
        Ok(Self::new(name, nuc_seq))
    }
    
    /// Convert an I/O Sequence to a ReferenceSequence
    pub fn from_io_sequence(seq: &IoSequence) -> Result<Self> {
        Ok(Self {
            sequence: seq.seq.clone(),
            name: seq.id.clone(),
            length: seq.len(),
        })
    }

    /// Get the length of the reference sequence
    pub fn len(&self) -> usize {
        self.length
    }

    /// Check if the reference sequence is empty
    pub fn is_empty(&self) -> bool {
        self.length == 0
    }

    /// Get the name of the reference sequence
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Get the base at a given position
    pub fn base_at(&self, pos: usize) -> Option<Nucleotide> {
        if pos < self.length {
            Some(self.sequence.get(pos))
        } else {
            None
        }
    }
    
    /// Get a slice of the sequence
    pub fn slice(&self, range: Range<usize>) -> Option<NucSeq> {
        if range.start <= range.end && range.end <= self.length {
            Some(self.sequence.subseq(range.start, range.end))
        } else {
            None
        }
    }
    
    /// Get a subsequence as a new ReferenceSequence
    pub fn subsequence(&self, range: Range<usize>) -> Option<Self> {
        let start = range.start;
        let end = range.end;
        self.slice(range).map(|seq| {
            Self::new(&format!("{}_{}_{}",  self.name, start, end), seq)
        })
    }

    /// Get the full sequence
    pub fn sequence(&self) -> &NucSeq {
        &self.sequence
    }
}

impl fmt::Display for ReferenceSequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, ">{}\n{}", self.name, self.sequence)
    }
}

/// Collection of reference sequences
pub type ReferenceSequences = HashMap<String, ReferenceSequence>; 