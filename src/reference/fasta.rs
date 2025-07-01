use std::path::Path;
use std::collections::HashMap;
use thiserror::Error;
use anyhow::{Result, Context};

use super::sequence::{ReferenceSequence, ReferenceSequences};
use crate::io::fasta as fasta_io;
// use crate::io::sequence::Sequence as IoSequence; // Unused for now

#[derive(Error, Debug)]
pub enum FastaError {
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),
    
    #[error("Invalid FASTA format: {0}")]
    FormatError(String),
    
    #[error("Reference not found: {0}")]
    NotFound(String),
}

/// Represents a reference genome from a FASTA file
#[derive(Debug, Clone)]
pub struct FastaReference {
    /// The sequences in the reference
    sequences: ReferenceSequences,
    /// The path to the FASTA file
    path: String,
}

impl FastaReference {
    /// Create a new FastaReference from a file path
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path_str = path.as_ref().to_string_lossy().to_string();
        
        // Use our new seq_io-based FASTA reader to load sequences
        let io_sequences = fasta_io::load_reference(&path)
            .with_context(|| format!("Failed to load reference from FASTA file: {}", path_str))?;
        
        // Convert IoSequence to ReferenceSequence
        let mut sequences = HashMap::new();
        
        for io_seq in io_sequences {
            let ref_seq = ReferenceSequence::from_io_sequence(&io_seq)?;
            sequences.insert(ref_seq.name().to_string(), ref_seq);
        }
        
        Ok(Self {
            sequences,
            path: path_str,
        })
    }
    
    /// Get all sequences in the reference
    pub fn sequences(&self) -> &ReferenceSequences {
        &self.sequences
    }
    
    /// Get a specific sequence by name
    pub fn get_sequence(&self, name: &str) -> Option<&ReferenceSequence> {
        self.sequences.get(name)
    }
    
    /// Get the path to the FASTA file
    pub fn path(&self) -> &str {
        &self.path
    }
    
    /// Get the number of sequences in the reference
    pub fn len(&self) -> usize {
        self.sequences.len()
    }
    
    /// Check if the reference has any sequences
    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }
} 