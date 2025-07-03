use std::path::{Path, PathBuf};
use std::fs::File;
use std::io::BufWriter;
use anyhow::Result;
use serde::{Serialize, Deserialize};
use memmap2::Mmap;

use super::sequence::ReferenceSequence;

/// Configuration for a hash table
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HashtableConfig {
    /// K-mer size used for the hash table
    pub kmer_size: usize,
    /// Hash table type (normal, anchored, etc.)
    pub hash_type: HashTableType,
    /// Reference file path that was used to build the hash table
    pub reference_path: PathBuf,
    /// Output directory where hash table is stored
    pub output_dir: PathBuf,
    /// Number of threads used to build the hash table
    pub threads: usize,
    /// Version of the hash table format
    pub version: String,
}

/// Type of hash table
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum HashTableType {
    /// Normal hash table (default)
    Normal,
    /// Anchored hash table (for RNA-seq)
    Anchored,
    /// Methylation hash table (C to T conversion)
    MethylCtoT,
    /// Methylation hash table (G to A conversion)
    MethylGtoA,
}

impl HashtableConfig {
    /// Create a new HashtableConfig
    pub fn new(
        kmer_size: usize,
        hash_type: HashTableType,
        reference_path: PathBuf,
        output_dir: PathBuf,
        threads: usize,
    ) -> Self {
        Self {
            kmer_size,
            hash_type,
            reference_path,
            output_dir,
            threads,
            version: env!("CARGO_PKG_VERSION").to_string(),
        }
    }

    /// Save the configuration to a file
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let file = File::create(path)?;
        let writer = BufWriter::new(file);
        serde_json::to_writer_pretty(writer, self)?;
        Ok(())
    }

    /// Load the configuration from a file (DRAGMAP format)
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        use std::io::{BufRead, BufReader};
        use anyhow::anyhow;
        
        let file = File::open(&path)?;
        let reader = BufReader::new(file);
        
        let mut kmer_size = 21; // Default
        let mut reference_path = PathBuf::new();
        let mut output_dir = PathBuf::new();
        let mut threads = 1;
        
        for line in reader.lines() {
            let line = line?;
            let line = line.trim();
            
            // Skip comments and empty lines
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            
            // Parse key = value format
            if let Some((key, value)) = line.split_once('=') {
                let key = key.trim();
                let value = value.trim().trim_matches('\'').trim_matches('"');
                
                match key {
                    "pri_seed_bases" => kmer_size = value.parse()?,
                    "reference_source" => reference_path = PathBuf::from(value),
                    "hash_table" => {
                        // Extract output directory from hash table path
                        if let Some(parent) = PathBuf::from(value).parent() {
                            output_dir = parent.to_path_buf();
                        }
                    },
                    "num_threads" => threads = value.parse()?,
                    _ => {} // Ignore other keys
                }
            }
        }
        
        if reference_path.as_os_str().is_empty() {
            return Err(anyhow!("Missing reference_source in config file"));
        }
        
        if output_dir.as_os_str().is_empty() {
            output_dir = path.as_ref().parent().unwrap_or(Path::new(".")).to_path_buf();
        }
        
        Ok(Self::new(
            kmer_size,
            HashTableType::Normal,
            reference_path,
            output_dir,
            threads,
        ))
    }
}

/// The hashtable used for mapping
#[derive(Debug)]
pub struct Hashtable {
    /// The configuration for the hash table
    config: HashtableConfig,
    /// The actual hash table data
    data: HashtableData,
    /// Extended table data (optional)
    extend_data: Option<HashtableData>,
}

/// Enum representing the hash table data storage
#[derive(Debug)]
pub enum HashtableData {
    /// Memory-mapped data
    MemoryMapped(Mmap),
    /// In-memory data
    InMemory(Vec<u64>),
}

impl Hashtable {
    /// Create a new hash table from config and data
    pub fn new(
        config: HashtableConfig, 
        data: HashtableData, 
        extend_data: Option<HashtableData>
    ) -> Self {
        Self {
            config,
            data,
            extend_data,
        }
    }

    /// Build a new hash table from a reference sequence
    pub fn build(
        _reference: &ReferenceSequence,
        _kmer_size: usize,
        _hash_type: HashTableType,
        _threads: usize,
        _output_dir: PathBuf,
    ) -> Result<Self> {
        // This is a placeholder implementation
        // TODO: Implement the actual hash table generation algorithm
        unimplemented!("Hashtable generation not yet implemented");
    }

    /// Get a reference to the hash table configuration
    pub fn config(&self) -> &HashtableConfig {
        &self.config
    }

    /// Map a k-mer to the hash table and find potential reference positions
    pub fn map_kmer(&self, _kmer: &[u8]) -> Result<Vec<u64>> {
        // This is a placeholder implementation
        // TODO: Implement k-mer mapping
        unimplemented!("K-mer mapping not yet implemented");
    }

    /// Get the size of the hash table in bytes
    pub fn size_bytes(&self) -> usize {
        match &self.data {
            HashtableData::MemoryMapped(mmap) => mmap.len(),
            HashtableData::InMemory(vec) => vec.len() * std::mem::size_of::<u64>(),
        }
    }
} 