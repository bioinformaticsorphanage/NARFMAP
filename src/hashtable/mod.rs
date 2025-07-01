use std::path::{Path, PathBuf};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use anyhow::{Result, anyhow};
use log::{info, debug, warn};
use rayon::prelude::*;

use crate::config::{Config, HashTableConfig};
use crate::reference::sequence::{ReferenceSequence, NucSeq, Nucleotide};
use crate::reference::hashtable::{Hashtable as RefHashtable, HashtableConfig as RefHashtableConfig, HashTableType, HashtableData};

/// CRC polynomial constants for k-mer hashing
/// These match the CRC polynomials used in the original DRAGMAP implementation
const CRC_POLYNOMIALS: [u64; 4] = [
    0x1021, // CRC-16-CCITT
    0x8005, // CRC-16-IBM
    0x8408, // CRC-16-KERMIT  
    0x8810, // CRC-16-DNP
];

/// Hash table builder for generating hash tables from reference genomes
pub struct HashTableBuilder {
    config: HashTableConfig,
}

impl HashTableBuilder {
    /// Create a new hash table builder with the given configuration
    pub fn new(config: HashTableConfig) -> Self {
        Self { config }
    }

    /// Build a hash table from a reference FASTA file
    pub fn build_from_fasta<P: AsRef<Path>>(&self, fasta_path: P, output_dir: P) -> Result<RefHashtable> {
        let fasta_path = fasta_path.as_ref();
        let output_dir = output_dir.as_ref();

        info!("Building hash table from reference: {}", fasta_path.display());
        info!("Output directory: {}", output_dir.display());
        info!("K-mer size: {}", self.config.seed_len);
        info!("Using {} threads", self.config.num_threads);

        // Create output directory if it doesn't exist
        std::fs::create_dir_all(output_dir)?;

        // Load reference sequences
        let reference_sequences = crate::io::fasta::load_reference(fasta_path)?;
        info!("Loaded {} reference sequences", reference_sequences.len());

        // Build hash table for each sequence
        let mut hash_table_entries: HashMap<u64, Vec<u64>> = HashMap::new();
        let mut total_kmers = 0;

        for (seq_idx, sequence) in reference_sequences.iter().enumerate() {
            info!("Processing sequence {}: {} ({} bp)", 
                  seq_idx + 1, sequence.id, sequence.len());

            let kmers = self.extract_kmers(sequence, seq_idx as u64)?;
            total_kmers += kmers.len();

            // Add k-mers to hash table
            for (kmer_hash, positions) in kmers {
                hash_table_entries.entry(kmer_hash)
                    .or_insert_with(Vec::new)
                    .extend(positions);
            }
        }

        info!("Extracted {} total k-mers", total_kmers);
        info!("Hash table has {} unique k-mers", hash_table_entries.len());

        // Filter k-mers by frequency
        let filtered_entries = self.filter_by_frequency(hash_table_entries)?;
        info!("After frequency filtering: {} k-mers", filtered_entries.len());

        // Serialize hash table to disk
        let hash_table = self.serialize_hash_table(filtered_entries, fasta_path, output_dir)?;

        info!("Hash table construction complete");
        Ok(hash_table)
    }

    /// Extract k-mers from a reference sequence
    fn extract_kmers(&self, sequence: &crate::io::sequence::Sequence, seq_idx: u64) -> Result<HashMap<u64, Vec<u64>>> {
        let mut kmers: HashMap<u64, Vec<u64>> = HashMap::new();
        let k = self.config.seed_len;

        if sequence.len() < k {
            warn!("Sequence {} is shorter than k-mer size {}, skipping", sequence.id, k);
            return Ok(kmers);
        }

        // Extract k-mers in parallel chunks
        let chunk_size = 10000;
        let positions: Vec<usize> = (0..=sequence.len() - k).collect();
        
        let chunk_results: Vec<HashMap<u64, Vec<u64>>> = positions
            .par_chunks(chunk_size)
            .map(|chunk| {
                let mut chunk_kmers: HashMap<u64, Vec<u64>> = HashMap::new();
                
                for &pos in chunk {
                    if let Ok(kmer_seq) = NucSeq::from_str(&sequence.seq.to_string()[pos..pos + k]) {
                        if let Ok(kmer_hash) = kmer_seq.to_2bit_kmer(k) {
                            // Apply CRC hashing
                            let crc_hash = self.apply_crc_hash(kmer_hash, self.config.crc_primary);
                            
                            // Encode position with sequence index
                            let encoded_pos = (seq_idx << 32) | (pos as u64);
                            
                            chunk_kmers.entry(crc_hash)
                                .or_insert_with(Vec::new)
                                .push(encoded_pos);
                        }
                    }
                }
                
                chunk_kmers
            })
            .collect();

        // Merge results from all chunks
        for chunk_kmers in chunk_results {
            for (kmer_hash, positions) in chunk_kmers {
                kmers.entry(kmer_hash)
                    .or_insert_with(Vec::new)
                    .extend(positions);
            }
        }

        Ok(kmers)
    }

    /// Apply CRC hashing to a k-mer
    fn apply_crc_hash(&self, kmer: u64, crc_index: u32) -> u64 {
        if crc_index as usize >= CRC_POLYNOMIALS.len() {
            return kmer; // Fall back to identity if invalid CRC index
        }

        let polynomial = CRC_POLYNOMIALS[crc_index as usize];
        let mut hash = kmer;

        // Simple CRC-like hash transformation
        for _ in 0..16 {
            if hash & 0x8000000000000000 != 0 {
                hash = (hash << 1) ^ polynomial;
            } else {
                hash <<= 1;
            }
        }

        hash
    }

    /// Filter k-mers by frequency to reduce hash table size
    fn filter_by_frequency(&self, mut entries: HashMap<u64, Vec<u64>>) -> Result<HashMap<u64, Vec<u64>>> {
        let max_freq = self.config.max_seed_freq as usize;
        let target_freq = self.config.target_seed_freq as usize;

        debug!("Filtering k-mers with max frequency: {}", max_freq);

        // Remove k-mers that occur too frequently
        entries.retain(|_, positions| positions.len() <= max_freq);

        // For k-mers above target frequency, apply soft capping
        let soft_cap = self.config.soft_seed_freq_cap as usize;
        for (_, positions) in entries.iter_mut() {
            if positions.len() > soft_cap {
                // Keep only a subset of positions for high-frequency k-mers
                positions.truncate(target_freq);
            }
        }

        Ok(entries)
    }

    /// Serialize hash table to disk
    fn serialize_hash_table(
        &self,
        entries: HashMap<u64, Vec<u64>>,
        reference_path: &Path,
        output_dir: &Path,
    ) -> Result<RefHashtable> {
        // Create configuration for the reference hashtable
        let ref_config = RefHashtableConfig::new(
            self.config.seed_len,
            HashTableType::Normal,
            reference_path.to_path_buf(),
            output_dir.to_path_buf(),
            self.config.num_threads,
        );

        // Convert HashMap to Vec<u64> for storage
        let mut hash_data = Vec::new();
        
        // Simple serialization format: [kmer_hash, position_count, position1, position2, ...]
        for (kmer_hash, positions) in entries.iter() {
            hash_data.push(*kmer_hash);
            hash_data.push(positions.len() as u64);
            hash_data.extend(positions.iter());
        }

        info!("Serialized {} hash table entries into {} u64 values", 
              entries.len(), hash_data.len());

        // Create hash table data
        let data = HashtableData::InMemory(hash_data);

        // Save configuration to file
        let config_path = output_dir.join("hash_table.cfg");
        ref_config.save(&config_path)?;
        info!("Saved hash table configuration to: {}", config_path.display());

        // Create and return the hash table
        Ok(RefHashtable::new(ref_config, data, None))
    }
}

/// Hash table query interface for alignment
pub struct HashTableQuery {
    hashtable: RefHashtable,
    config: HashTableConfig,
}

impl HashTableQuery {
    /// Create a new hash table query interface
    pub fn new(hashtable: RefHashtable, config: HashTableConfig) -> Self {
        Self { hashtable, config }
    }

    /// Load a hash table from a directory
    pub fn load_from_dir<P: AsRef<Path>>(dir_path: P, config: HashTableConfig) -> Result<Self> {
        let dir_path = dir_path.as_ref();
        let config_path = dir_path.join("hash_table.cfg");

        if !config_path.exists() {
            return Err(anyhow!("Hash table configuration not found: {}", config_path.display()));
        }

        // Load configuration
        let ref_config = RefHashtableConfig::load(&config_path)?;
        
        // For now, create empty hash table (actual loading will be implemented later)
        let data = HashtableData::InMemory(Vec::new());
        let hashtable = RefHashtable::new(ref_config, data, None);

        Ok(Self::new(hashtable, config))
    }

    /// Query the hash table for potential positions of a k-mer
    pub fn query_kmer(&self, kmer: &[u8]) -> Result<Vec<u64>> {
        if kmer.len() != self.config.seed_len {
            return Err(anyhow!("K-mer length {} does not match expected length {}", 
                              kmer.len(), self.config.seed_len));
        }

        // Convert k-mer to hash value
        let kmer_seq = NucSeq::from_bytes(kmer);
        let kmer_hash = kmer_seq.to_2bit_kmer(self.config.seed_len)?;
        
        // Apply same CRC transformation as during building
        let crc_hash = self.apply_crc_hash(kmer_hash, self.config.crc_primary);

        // Query the hash table (placeholder implementation)
        // TODO: Implement actual hash table lookup
        Ok(Vec::new())
    }

    /// Apply CRC hashing to a k-mer (same as in builder)
    fn apply_crc_hash(&self, kmer: u64, crc_index: u32) -> u64 {
        if crc_index as usize >= CRC_POLYNOMIALS.len() {
            return kmer;
        }

        let polynomial = CRC_POLYNOMIALS[crc_index as usize];
        let mut hash = kmer;

        for _ in 0..16 {
            if hash & 0x8000000000000000 != 0 {
                hash = (hash << 1) ^ polynomial;
            } else {
                hash <<= 1;
            }
        }

        hash
    }

    /// Get hash table statistics
    pub fn get_stats(&self) -> HashMap<String, u64> {
        let mut stats = HashMap::new();
        stats.insert("size_bytes".to_string(), self.hashtable.size_bytes() as u64);
        stats.insert("kmer_size".to_string(), self.config.seed_len as u64);
        stats.insert("max_seed_freq".to_string(), self.config.max_seed_freq as u64);
        stats
    }
}

/// Legacy Hashtable struct for compatibility
/// This is kept for backward compatibility with existing code
#[deprecated(note = "Use HashTableBuilder and HashTableQuery instead")]
pub struct Hashtable {
    query: Option<HashTableQuery>,
}

#[allow(deprecated)]
impl Hashtable {
    pub fn new() -> Self {
        Self { query: None }
    }

    /// Create from a config (for compatibility)
    pub fn from_config(config: HashTableConfig) -> Self {
        Self { query: None }
    }
} 