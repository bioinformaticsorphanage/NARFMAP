pub mod sequence;
pub mod fasta;
pub mod hashtable;

use std::path::Path;
use anyhow::Result;

pub use sequence::ReferenceSequence;
pub use hashtable::{HashtableConfig, Hashtable};
pub use fasta::FastaReference;

/// Loads a reference genome from a FASTA file
pub fn load_reference<P: AsRef<Path>>(path: P) -> Result<FastaReference> {
    FastaReference::from_path(path)
}

/// Interface for the reference directory containing the hash table
#[derive(Debug)]
pub struct ReferenceDir {
    reference_sequence: ReferenceSequence,
    hashtable_config: HashtableConfig,
    hashtable: Hashtable,
}

impl ReferenceDir {
    /// Create a new ReferenceDir from a directory path
    pub fn new<P: AsRef<Path>>(path: P, mmap_reference: bool, load_reference: bool) -> Result<Self> {
        // This is a placeholder implementation that will be filled in later
        unimplemented!("ReferenceDir::new not yet implemented")
    }

    /// Get the reference sequence
    pub fn get_reference_sequence(&self) -> &ReferenceSequence {
        &self.reference_sequence
    }

    /// Get the hashtable configuration
    pub fn get_hashtable_config(&self) -> &HashtableConfig {
        &self.hashtable_config
    }
    
    /// Get the hashtable
    pub fn get_hashtable(&self) -> &Hashtable {
        &self.hashtable
    }
} 