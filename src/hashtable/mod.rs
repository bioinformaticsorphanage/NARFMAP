use std::path::{Path, PathBuf};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use anyhow::{Result, anyhow};
use log::{info, debug, warn};
use rayon::prelude::*;

use crate::config::{Config, HashTableConfig};
use crate::reference::hashtable::{Hashtable as RefHashtable, HashtableConfig as RefHashtableConfig, HashTableType, HashtableData};

mod crc_hash;
mod hash_record;
mod bucket;

pub use crc_hash::{CrcPolynomial, CrcHasher, KmerHasher};
pub use hash_record::{HashRecord, RecordType};
pub use bucket::{Bucket, HashtableTraits};

/// Hash table builder for generating hash tables from reference genomes
pub struct HashTableBuilder {
    config: HashTableConfig,
    kmer_hasher: KmerHasher,
}

impl HashTableBuilder {
    /// Create a new hash table builder with the given configuration
    pub fn new(config: HashTableConfig) -> Result<Self> {
        let kmer_hasher = KmerHasher::with_dragmap_defaults(config.seed_len)?;
        Ok(Self { config, kmer_hasher })
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

    /// Extract k-mers from a reference sequence using efficient 2-bit encoding
    /// This matches DRAGMAP's approach for k-mer extraction and position encoding
    fn extract_kmers(&self, sequence: &crate::io::sequence::Sequence, seq_idx: u64) -> Result<HashMap<u64, Vec<u64>>> {
        let mut kmers: HashMap<u64, Vec<u64>> = HashMap::new();
        let k = self.config.seed_len;

        if sequence.len() < k {
            warn!("Sequence {} is shorter than k-mer size {}, skipping", sequence.id, k);
            return Ok(kmers);
        }

        // Use DRAGMAP-style seed interval for sampling
        let interval = self.config.ref_seed_interval.max(1.0) as usize;
        
        // Extract k-mers using efficient sliding window approach
        let chunk_size = 10000;
        let total_positions = sequence.len() - k + 1;
        let positions: Vec<usize> = (0..total_positions).step_by(interval).collect();
        
        debug!("Extracting k-mers from sequence {} ({} bp) with interval {}, {} positions", 
               sequence.id, sequence.len(), interval, positions.len());
        
        let chunk_results: Vec<HashMap<u64, Vec<u64>>> = positions
            .par_chunks(chunk_size)
            .map(|chunk| {
                let mut chunk_kmers: HashMap<u64, Vec<u64>> = HashMap::new();
                
                for &pos in chunk {
                    // Extract k-mer using direct 2-bit access for efficiency
                    if let Some(kmer_subseq) = sequence.substring(pos, pos + k) {
                        // Convert to 2-bit packed representation
                        if let Ok(kmer_2bit) = self.extract_kmer_2bit(&kmer_subseq) {
                            // Hash the k-mer using CRC
                            if let Ok(crc_hash) = self.kmer_hasher.hash_kmer(kmer_2bit) {
                                // Encode position with sequence index (DRAGMAP format)
                                // Upper 32 bits: sequence index, lower 32 bits: position
                                let encoded_pos = (seq_idx << 32) | (pos as u64);
                                
                                chunk_kmers.entry(crc_hash)
                                    .or_insert_with(Vec::new)
                                    .push(encoded_pos);
                            }
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

        debug!("Extracted {} unique k-mers from sequence {}", kmers.len(), sequence.id);
        Ok(kmers)
    }

    /// Extract k-mer as 2-bit packed representation efficiently
    /// This avoids string conversion and matches DRAGMAP's 2-bit encoding
    fn extract_kmer_2bit(&self, kmer_seq: &crate::reference::sequence::NucSeq) -> Result<u64> {
        if kmer_seq.len() > 32 {
            return Err(anyhow!("K-mer too long for u64 encoding: {} bases", kmer_seq.len()));
        }
        
        let mut kmer_2bit = 0u64;
        for i in 0..kmer_seq.len() {
            let base_2bit = match kmer_seq.get(i) {
                crate::reference::sequence::Nucleotide::A => 0u64,
                crate::reference::sequence::Nucleotide::C => 1u64, 
                crate::reference::sequence::Nucleotide::G => 2u64,
                crate::reference::sequence::Nucleotide::T => 3u64,
                crate::reference::sequence::Nucleotide::N => 2u64, // N treated as G (DRAGMAP compatible)
            };
            
            // Pack in little-endian format (2 bits per base)
            kmer_2bit |= base_2bit << (i * 2);
        }
        
        Ok(kmer_2bit)
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

    /// Build bucket-based hash table from k-mer entries
    fn build_bucket_table(&self, entries: HashMap<u64, Vec<u64>>) -> Result<Vec<Bucket>> {
        // Calculate hash table size based on CRC bits
        let table_size_bits = self.config.crc_primary as usize;
        let squeeze_factor = 1u64; // Start with no squeeze
        
        // Ensure minimum table size to avoid underflow
        let effective_bits = table_size_bits.max(10); // Minimum 10 bits
        let num_buckets = (1u64 << (effective_bits - 6)) as usize; // 2^(bits-6) buckets (64 bytes each)
        
        info!("Building hash table with {} buckets ({} MB)", 
              num_buckets, (num_buckets * 64) / (1024 * 1024));
        
        // Initialize buckets
        let mut buckets: Vec<Bucket> = vec![Bucket::new(); num_buckets];
        let mut collision_count = 0u64;
        let mut chain_count = 0u64;
        
        // Process each k-mer and its positions
        for (kmer_hash, positions) in entries.iter() {
            // Get virtual address and bucket index
            let virtual_addr = self.kmer_hasher.get_address_from_hash(*kmer_hash, squeeze_factor);
            let bucket_idx = (virtual_addr >> HashtableTraits::HASH_BUCKET_BYTES_LOG2) as usize % num_buckets;
            
            // Extract hash bits for the record (bits 35-57)
            let hash_bits = ((*kmer_hash >> 35) & 0x7FFFFF) as u32;
            
            // Get thread ID from virtual address
            let thread_id = self.kmer_hasher.get_thread_id_from_address(virtual_addr);
            
            // Handle high-frequency k-mers
            if positions.len() > self.config.max_seed_freq as usize {
                // Create HIFREQ record
                let hifreq_record = HashRecord::hifreq(
                    thread_id,
                    hash_bits,
                    false, // not extended
                    false, // not last
                    false, // no random sample
                    false, // not alt
                    positions.len() as u32
                );
                
                // Try to insert in bucket
                if !self.insert_record_in_bucket(&mut buckets[bucket_idx], hifreq_record) {
                    collision_count += 1;
                    // TODO: Implement linear probing or chaining
                }
                continue;
            }
            
            // Create HIT records for each position
            for (i, &encoded_pos) in positions.iter().enumerate() {
                let seq_idx = (encoded_pos >> 32) as u32;
                let position = (encoded_pos & 0xFFFFFFFF) as u32;
                let is_last = i == positions.len() - 1;
                
                let hit_record = HashRecord::hit(
                    thread_id,
                    hash_bits,
                    false, // not extended (primary seeds)
                    is_last,
                    false, // TODO: detect reverse complement
                    position
                );
                
                // Try to insert in bucket
                if !self.insert_record_in_bucket(&mut buckets[bucket_idx], hit_record) {
                    collision_count += 1;
                    // TODO: Implement linear probing for overflow
                }
            }
        }
        
        info!("Hash table built: {} collisions, {} chains", collision_count, chain_count);
        Ok(buckets)
    }
    
    /// Insert a record into a bucket, returns true if successful
    fn insert_record_in_bucket(&self, bucket: &mut Bucket, record: HashRecord) -> bool {
        if let Some(slot) = bucket.find_empty_slot() {
            bucket.set(slot, record).unwrap();
            true
        } else {
            false
        }
    }

    /// Serialize hash table to disk
    fn serialize_hash_table(
        &self,
        entries: HashMap<u64, Vec<u64>>,
        reference_path: &Path,
        output_dir: &Path,
    ) -> Result<RefHashtable> {
        // Build bucket-based hash table
        let buckets = self.build_bucket_table(entries)?;
        
        // Create binary file for hash table
        let hash_table_path = output_dir.join("hash_table.bin");
        let mut file = std::fs::File::create(&hash_table_path)?;
        
        // Write buckets to file
        use std::io::Write;
        for bucket in buckets.iter() {
            file.write_all(&bucket.to_bytes())?;
        }
        
        info!("Wrote hash table to: {}", hash_table_path.display());
        
        // Create configuration for the reference hashtable
        let ref_config = RefHashtableConfig::new(
            self.config.seed_len,
            HashTableType::Normal,
            reference_path.to_path_buf(),
            output_dir.to_path_buf(),
            self.config.num_threads,
        );

        // Save configuration to file
        let config_path = output_dir.join("hash_table.cfg");
        ref_config.save(&config_path)?;
        info!("Saved hash table configuration to: {}", config_path.display());

        // Create hash table data (using empty in-memory for now)
        let data = HashtableData::InMemory(Vec::new());

        // Create and return the hash table
        Ok(RefHashtable::new(ref_config, data, None))
    }
}

/// Hash table query interface for alignment
pub struct HashTableQuery {
    hashtable: RefHashtable,
    config: HashTableConfig,
    kmer_hasher: KmerHasher,
    buckets: Vec<Bucket>,
    num_buckets: usize,
}

impl HashTableQuery {
    /// Create a new hash table query interface with loaded buckets
    pub fn new(hashtable: RefHashtable, config: HashTableConfig, buckets: Vec<Bucket>) -> Result<Self> {
        let kmer_hasher = KmerHasher::with_dragmap_defaults(config.seed_len)?;
        let num_buckets = buckets.len();
        Ok(Self { hashtable, config, kmer_hasher, buckets, num_buckets })
    }

    /// Load a hash table from a directory
    pub fn load_from_dir<P: AsRef<Path>>(dir_path: P, config: HashTableConfig) -> Result<Self> {
        let dir_path = dir_path.as_ref();
        let config_path = dir_path.join("hash_table.cfg");
        let bin_path = dir_path.join("hash_table.bin");

        if !config_path.exists() {
            return Err(anyhow!("Hash table configuration not found: {}", config_path.display()));
        }

        if !bin_path.exists() {
            return Err(anyhow!("Hash table binary not found: {}", bin_path.display()));
        }

        // Load configuration
        let ref_config = RefHashtableConfig::load(&config_path)?;
        info!("Loaded hash table config: k-mer size {}", ref_config.kmer_size);
        
        // Load binary hash table data
        let buckets = Self::load_buckets_from_file(&bin_path, &config)?;
        info!("Loaded {} buckets from {}", buckets.len(), bin_path.display());
        
        // Create hash table with loaded data
        let data = HashtableData::InMemory(Vec::new()); // We use buckets instead
        let hashtable = RefHashtable::new(ref_config, data, None);

        Self::new(hashtable, config, buckets)
    }

    /// Load buckets from the binary hash table file
    fn load_buckets_from_file<P: AsRef<Path>>(bin_path: P, config: &HashTableConfig) -> Result<Vec<Bucket>> {
        use std::fs::File;
        use std::io::Read;
        
        let mut file = File::open(bin_path.as_ref())?;
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)?;
        
        debug!("Read {} bytes from hash table file", buffer.len());
        
        // Each bucket is 64 bytes (8 records * 8 bytes each)
        const BUCKET_SIZE: usize = 64;
        if buffer.len() % BUCKET_SIZE != 0 {
            return Err(anyhow!("Invalid hash table file: size {} is not divisible by bucket size {}", 
                              buffer.len(), BUCKET_SIZE));
        }
        
        let num_buckets = buffer.len() / BUCKET_SIZE;
        let mut buckets = Vec::with_capacity(num_buckets);
        
        // Parse buckets from binary data
        for i in 0..num_buckets {
            let start = i * BUCKET_SIZE;
            let end = start + BUCKET_SIZE;
            let bucket_bytes: [u8; 64] = buffer[start..end].try_into()
                .map_err(|_| anyhow!("Failed to extract bucket {} from binary data", i))?;
            
            buckets.push(Bucket::from_bytes(&bucket_bytes));
        }
        
        info!("Successfully loaded {} buckets", num_buckets);
        Ok(buckets)
    }

    /// Query the hash table for potential positions of a k-mer
    pub fn query_kmer(&self, kmer: &[u8]) -> Result<Vec<u64>> {
        if kmer.len() != self.config.seed_len {
            return Err(anyhow!("K-mer length {} does not match expected length {}", 
                              kmer.len(), self.config.seed_len));
        }

        // Convert k-mer to 2-bit representation and hash with CRC
        let kmer_2bit = self.kmer_hasher.sequence_to_2bit(kmer)?;
        let crc_hash = self.kmer_hasher.hash_kmer(kmer_2bit)?;

        debug!("Querying k-mer hash {:016x} in {} buckets", crc_hash, self.num_buckets);

        // Calculate bucket index using the same logic as in build_bucket_table
        let squeeze_factor = 1u64;
        let virtual_addr = self.kmer_hasher.get_address_from_hash(crc_hash, squeeze_factor);
        let bucket_idx = (virtual_addr >> HashtableTraits::HASH_BUCKET_BYTES_LOG2) as usize % self.num_buckets;
        
        // Extract the hash bits we're looking for (bits 35-57, same as in builder)
        let target_hash_bits = ((crc_hash >> 35) & 0x7FFFFF) as u32;
        
        debug!("Looking in bucket {} for hash bits {:06x}", bucket_idx, target_hash_bits);
        
        let mut positions = Vec::new();
        
        if bucket_idx < self.buckets.len() {
            let bucket = &self.buckets[bucket_idx];
            
            // Search through all records in the bucket
            for i in 0..8 { // 8 records per bucket
                if let Some(record) = bucket.get(i) {
                    // Check if this record matches our target hash
                    if record.hash_bits() == target_hash_bits {
                        match record.record_type() {
                            RecordType::Hit => {
                                if let Some(ref_pos) = record.reference_position() {
                                    // Reconstruct the encoded position (sequence_id in upper 32 bits, position in lower 32)
                                    // For now, assume sequence_id = 0 since we don't store it in the record
                                    let encoded_pos = (0u64 << 32) | (ref_pos as u64);
                                    positions.push(encoded_pos);
                                    
                                    debug!("Found HIT: position {}", ref_pos);
                                    
                                    // If this is the last record in a chain, stop
                                    if record.is_last_in_thread() {
                                        break;
                                    }
                                }
                            }
                            RecordType::HiFreq => {
                                debug!("Found HIFREQ record with frequency {}", record.frequency().unwrap_or(0));
                                // For high-frequency k-mers, we might skip them or handle specially
                                // For now, just return empty to avoid overwhelming results
                                break;
                            }
                            _ => {
                                debug!("Found other record type: {:?}", record.record_type());
                            }
                        }
                    }
                }
            }
        }
        
        debug!("Found {} positions for k-mer", positions.len());
        Ok(positions)
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
    pub fn from_config(_config: HashTableConfig) -> Self {
        Self { query: None }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    #[test]
    fn test_kmer_extraction_efficiency() {
        // Test the new efficient k-mer extraction
        let mut config = HashTableConfig::default();
        config.seed_len = 16;
        config.ref_seed_interval = 2.0; // Sample every 2 bases
        config.max_seed_freq = 1000;
        config.target_seed_freq = 100.0;
        config.soft_seed_freq_cap = 200.0;
        config.num_threads = 1;
        
        let builder = HashTableBuilder::new(config).unwrap();
        
        // Create a test sequence
        let test_seq = crate::io::sequence::Sequence::new(
            "test_seq".to_string(),
            "ACGTACGTACGTACGTACGTACGTACGTACGT", // 32 bases
            None
        ).unwrap();
        
        // Extract k-mers
        let kmers = builder.extract_kmers(&test_seq, 0).unwrap();
        
        // Should have extracted some k-mers (with interval=2, from 32-base seq with 16-mers)
        // Positions: 0, 2, 4, 6, 8, 10, 12, 14, 16 (9 positions)
        assert!(kmers.len() > 0);
        assert!(kmers.len() <= 9); // Max 9 due to sampling interval
        
        // Verify positions are encoded correctly (seq_idx=0 in upper 32 bits)
        for positions in kmers.values() {
            for &pos in positions {
                let seq_idx = pos >> 32;
                let position = pos & 0xFFFFFFFF;
                assert_eq!(seq_idx, 0);
                assert!(position < 32); // Position should be within sequence
            }
        }
    }
    
    #[test]
    fn test_2bit_kmer_extraction() {
        let mut config = HashTableConfig::default();
        config.seed_len = 4;
        config.ref_seed_interval = 1.0;
        config.max_seed_freq = 1000;
        config.target_seed_freq = 100.0;
        config.soft_seed_freq_cap = 200.0;
        config.num_threads = 1;
        
        let builder = HashTableBuilder::new(config).unwrap();
        
        // Test the direct 2-bit extraction method
        let test_nuc_seq = crate::reference::sequence::NucSeq::from_str("ACGT").unwrap();
        let kmer_2bit = builder.extract_kmer_2bit(&test_nuc_seq).unwrap();
        
        // ACGT = A(0) C(1) G(2) T(3) in 2-bit: 0b11100100 = 0xE4
        assert_eq!(kmer_2bit & 0xFF, 0b11100100);
    }
    
    #[test] 
    fn test_load_test_fasta() {
        // Test loading the actual tiny.fasta test file if it exists
        let test_fasta_path = Path::new("data/tiny/tiny-2x1Xrepeats.v8/tiny.fasta");
        
        if test_fasta_path.exists() {
            let sequences = crate::io::fasta::load_reference(test_fasta_path).unwrap();
            assert!(sequences.len() > 0);
            
            let first_seq = &sequences[0];
            assert!(first_seq.len() > 100); // Should have a decent length
            assert!(first_seq.id.contains("phiX174")); // Should be phiX174
            
            // Test k-mer extraction on real data
            let mut config = HashTableConfig::default();
            config.seed_len = 21;
            config.ref_seed_interval = 4.0; // Sample every 4 bases
            config.max_seed_freq = 1000;
            config.target_seed_freq = 100.0;
            config.soft_seed_freq_cap = 200.0;
            config.num_threads = 1;
            
            let builder = HashTableBuilder::new(config).unwrap();
            let kmers = builder.extract_kmers(first_seq, 0).unwrap();
            
            assert!(kmers.len() > 0);
            println!("Extracted {} unique k-mers from {} bp sequence", 
                     kmers.len(), first_seq.len());
        } else {
            println!("Test FASTA file not found, skipping integration test");
        }
    }
    
    #[test]
    fn test_bucket_table_construction() {
        let mut config = HashTableConfig::default();
        config.seed_len = 16;
        config.crc_primary = 16; // Small table for testing
        config.max_seed_freq = 10;
        
        let builder = HashTableBuilder::new(config).unwrap();
        
        // Create some test k-mer entries
        let mut entries = HashMap::new();
        
        // Add a normal frequency k-mer
        entries.insert(0x123456789ABCDEF0, vec![0x100, 0x200, 0x300]);
        
        // Add a high-frequency k-mer
        entries.insert(0xFEDCBA9876543210, vec![0x1000; 15]); // 15 positions > max_seed_freq
        
        // Build bucket table
        let buckets = builder.build_bucket_table(entries).unwrap();
        
        // Should have 2^(16-6) = 1024 buckets
        assert_eq!(buckets.len(), 1024);
        
        // Check that some buckets have records
        let non_empty_count = buckets.iter()
            .filter(|b| b.count_occupied() > 0)
            .count();
        assert!(non_empty_count > 0);
    }
    
    #[test]
    fn test_hash_record_in_bucket() {
        let mut bucket = Bucket::new();
        
        // Add different types of records
        let hit1 = HashRecord::hit(1, 0x12345, false, false, false, 1000);
        let hit2 = HashRecord::hit(2, 0x54321, true, true, true, 2000);
        let hifreq = HashRecord::hifreq(3, 0xABCDE, false, false, false, false, 100);
        
        bucket.set(0, hit1).unwrap();
        bucket.set(1, hit2).unwrap();
        bucket.set(2, hifreq).unwrap();
        
        // Verify records
        assert_eq!(bucket[0].record_type(), RecordType::Hit);
        assert_eq!(bucket[0].reference_position(), Some(1000));
        
        assert_eq!(bucket[1].record_type(), RecordType::Hit);
        assert!(bucket[1].is_extended());
        assert!(bucket[1].is_last_in_thread());
        assert!(bucket[1].is_reverse_complement());
        
        assert_eq!(bucket[2].record_type(), RecordType::HiFreq);
        assert_eq!(bucket[2].frequency(), Some(100));
    }
} 