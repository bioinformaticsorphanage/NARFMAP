use std::path::Path;
use std::collections::HashMap;
// use std::fs::File;
// use std::io::{BufWriter, Write}; // Unused for now
use anyhow::{Result, anyhow};
use log::{info, debug, warn};
use rayon::prelude::*;

use crate::config::HashTableConfig;
use crate::reference::hashtable::{Hashtable as RefHashtable, HashtableConfig as RefHashtableConfig, HashTableType, HashtableData};

mod crc_hash;
mod hash_record;
mod bucket;
mod reference_metadata;

pub use crc_hash::{CrcPolynomial, CrcHasher, KmerHasher};
pub use hash_record::{HashRecord, RecordType};
pub use bucket::{Bucket, HashtableTraits};
pub use reference_metadata::{ReferenceMetadata, SequenceInfo};

/// Hash table builder for generating hash tables from reference genomes
pub struct HashTableBuilder {
    config: HashTableConfig,
    kmer_hasher: KmerHasher,
    collision_stats: CollisionStats,
}

/// Statistics for hash table collision resolution
#[derive(Debug, Default)]
pub struct CollisionStats {
    pub total_insertions: u64,
    pub direct_insertions: u64,
    pub probing_insertions: u64,
    pub chaining_insertions: u64,
    pub failed_insertions: u64,
    pub total_probes: u64,
    pub max_probe_depth: u32,
    pub total_chains: u64,
    pub max_chain_length: u32,
}

impl HashTableBuilder {
    /// Create a new hash table builder with the given configuration
    pub fn new(config: HashTableConfig) -> Result<Self> {
        let kmer_hasher = KmerHasher::with_dragmap_defaults(config.seed_len)?;
        Ok(Self { 
            config, 
            kmer_hasher, 
            collision_stats: CollisionStats::default(),
        })
    }

    /// Build a hash table from a reference FASTA file
    pub fn build_from_fasta<P: AsRef<Path>>(&mut self, fasta_path: P, output_dir: P) -> Result<RefHashtable> {
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
        let hash_table = self.serialize_hash_table(
            filtered_entries, 
            &reference_sequences,
            fasta_path, 
            output_dir
        )?;

        info!("Hash table construction complete");
        self.print_collision_statistics();
        Ok(hash_table)
    }
    
    /// Print detailed collision resolution statistics
    fn print_collision_statistics(&self) {
        let stats = &self.collision_stats;
        info!("=== Hash Table Collision Statistics ===");
        info!("Total insertions: {}", stats.total_insertions);
        info!("Direct insertions: {} ({:.1}%)", 
              stats.direct_insertions, 
              stats.direct_insertions as f64 / stats.total_insertions as f64 * 100.0);
        info!("Probing insertions: {} ({:.1}%)", 
              stats.probing_insertions,
              stats.probing_insertions as f64 / stats.total_insertions as f64 * 100.0);
        info!("Chaining insertions: {} ({:.1}%)", 
              stats.chaining_insertions,
              stats.chaining_insertions as f64 / stats.total_insertions as f64 * 100.0);
        info!("Failed insertions: {} ({:.1}%)", 
              stats.failed_insertions,
              stats.failed_insertions as f64 / stats.total_insertions as f64 * 100.0);
        info!("Total probes: {}", stats.total_probes);
        info!("Max probe depth: {}", stats.max_probe_depth);
        info!("Total chains created: {}", stats.total_chains);
        
        if stats.total_insertions > 0 {
            let avg_probes = stats.total_probes as f64 / stats.total_insertions as f64;
            info!("Average probes per insertion: {:.2}", avg_probes);
        }
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
    fn build_bucket_table(&mut self, entries: HashMap<u64, Vec<u64>>) -> Result<Vec<Bucket>> {
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
        let chain_count = 0u64;
        
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
                
                // Try to insert in bucket with linear probing
                if !self.insert_record_with_probing(&mut buckets, bucket_idx, hifreq_record) {
                    collision_count += 1;
                    warn!("Failed to insert HIFREQ record after linear probing");
                }
                continue;
            }
            
            // Create HIT records for each position
            for (i, &encoded_pos) in positions.iter().enumerate() {
                let _seq_idx = (encoded_pos >> 32) as u32;
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
                
                // Try to insert in bucket with linear probing
                if !self.insert_record_with_probing(&mut buckets, bucket_idx, hit_record) {
                    collision_count += 1;
                    debug!("Collision for position {} in bucket {}", position, bucket_idx);
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
    
    /// Insert a record with DRAGMAP-compatible block-constrained probing and chaining
    /// This implements the sophisticated collision resolution used by DRAGMAP
    fn insert_record_with_probing(&mut self, buckets: &mut [Bucket], start_idx: usize, record: HashRecord) -> bool {
        self.collision_stats.total_insertions += 1;
        let num_buckets = buckets.len();
        
        // Calculate block constraints based on DRAGMAP's approach
        let buckets_per_block = self.get_buckets_per_block();
        let block_start = (start_idx / buckets_per_block) * buckets_per_block;
        let block_end = (block_start + buckets_per_block).min(num_buckets);
        
        debug!("Inserting record in block [{}, {}), start_idx: {}", block_start, block_end, start_idx);
        
        // First, try block-constrained linear probing (up to 8 buckets within block)
        let max_probes = 8.min(block_end - block_start);
        
        for probe in 0..max_probes {
            let bucket_idx = block_start + ((start_idx - block_start + probe) % (block_end - block_start));
            
            self.collision_stats.total_probes += 1;
            self.collision_stats.max_probe_depth = self.collision_stats.max_probe_depth.max(probe as u32 + 1);
            
            // Try to insert in current bucket
            if self.insert_record_in_bucket(&mut buckets[bucket_idx], record) {
                if probe == 0 {
                    self.collision_stats.direct_insertions += 1;
                } else {
                    self.collision_stats.probing_insertions += 1;
                    debug!("Inserted record in bucket {} after {} probes (block-constrained)", bucket_idx, probe);
                }
                return true;
            }
        }
        
        // If block-constrained probing fails, try chaining
        debug!("Block-constrained probing failed, attempting chaining");
        
        // Look for an available bucket to chain to (preferably in a different block)
        if let Some(chain_bucket_idx) = self.find_chain_bucket(buckets, block_start, block_end) {
            // Insert the record in the chain bucket
            if self.insert_record_in_bucket(&mut buckets[chain_bucket_idx], record) {
                debug!("Inserted record in chain bucket {}", chain_bucket_idx);
                
                // Create chain records in the original bucket
                if self.create_chain_records(buckets, start_idx, chain_bucket_idx, record) {
                    self.collision_stats.chaining_insertions += 1;
                    self.collision_stats.total_chains += 1;
                    return true;
                } else {
                    // If we can't create chain records, remove the record we just inserted
                    warn!("Failed to create chain records, cleaning up");
                    // TODO: Remove the record from chain_bucket_idx
                }
            }
        }
        
        debug!("All collision resolution attempts failed for bucket {}", start_idx);
        self.collision_stats.failed_insertions += 1;
        false
    }
    
    /// Calculate buckets per block based on CRC bit configuration
    /// This follows DRAGMAP's block organization for efficient probing
    fn get_buckets_per_block(&self) -> usize {
        // Use a reasonable block size based on the CRC primary bits
        // For good collision distribution, we want blocks to be large enough
        // but not so large that probing becomes inefficient
        let block_size_log2 = if self.config.crc_primary >= 16 {
            // For larger hash tables, use 256-bucket blocks (2^8 buckets = 2^14 bytes)
            8
        } else if self.config.crc_primary >= 12 {
            // For medium hash tables, use 64-bucket blocks (2^6 buckets = 2^12 bytes)
            6
        } else {
            // For small hash tables, use 16-bucket blocks (2^4 buckets = 2^10 bytes)
            4
        };
        
        let buckets_per_block_log2 = block_size_log2.max(3); // Minimum 8 buckets per block
        1 << buckets_per_block_log2
    }
    
    /// Find a suitable bucket for chaining outside the current block
    fn find_chain_bucket(&self, buckets: &[Bucket], block_start: usize, block_end: usize) -> Option<usize> {
        let num_buckets = buckets.len();
        let buckets_per_block = self.get_buckets_per_block();
        
        // Try to find a bucket in a different block with available space
        for block_offset in 1..=(num_buckets / buckets_per_block) {
            // Try both directions
            for direction in [1, -1] {
                let target_block_start = if direction > 0 {
                    (block_start + block_offset * buckets_per_block) % num_buckets
                } else {
                    (block_start + num_buckets - (block_offset * buckets_per_block)) % num_buckets
                };
                
                let target_block_end = (target_block_start + buckets_per_block).min(num_buckets);
                
                // Look for an available bucket in this block
                for bucket_idx in target_block_start..target_block_end {
                    if !buckets[bucket_idx].is_full() {
                        return Some(bucket_idx);
                    }
                }
            }
        }
        
        None
    }
    
    /// Create chain records to link the original bucket to the chain bucket
    /// This implements DRAGMAP's sophisticated chaining mechanism
    fn create_chain_records(&self, buckets: &mut [Bucket], original_bucket_idx: usize, chain_bucket_idx: usize, _original_record: HashRecord) -> bool {
        // Calculate chain pointer (relative bucket offset)
        let chain_pointer = if chain_bucket_idx > original_bucket_idx {
            chain_bucket_idx - original_bucket_idx
        } else {
            // Handle wrap-around
            buckets.len() + chain_bucket_idx - original_bucket_idx
        };
        
        // Ensure chain pointer fits in 18 bits
        if chain_pointer >= (1 << 18) {
            warn!("Chain pointer {} too large for 18-bit field", chain_pointer);
            return false;
        }
        
        // For simplicity, create a CHAIN_BEG_MASK record
        // In a full implementation, we would use filter masks and lists for efficiency
        let chain_record = HashRecord::chain_beg_mask(0xFFFFFFFF, chain_pointer as u32);
        
        // Try to insert the chain record in the original bucket
        if let Some(slot) = buckets[original_bucket_idx].find_empty_slot() {
            if buckets[original_bucket_idx].set(slot, chain_record).is_ok() {
                debug!("Created chain record in bucket {} pointing to bucket {}", original_bucket_idx, chain_bucket_idx);
                return true;
            }
        }
        
        warn!("No space for chain record in bucket {}", original_bucket_idx);
        false
    }

    /// Serialize hash table to disk
    fn serialize_hash_table(
        &mut self,
        entries: HashMap<u64, Vec<u64>>,
        reference_sequences: &[crate::io::sequence::Sequence],
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
        
        // Create reference metadata
        let ref_metadata = ReferenceMetadata::from_sequences(reference_sequences);
        
        // Save reference metadata
        let metadata_path = output_dir.join("reference_metadata.json");
        ref_metadata.save(&metadata_path)?;
        info!("Saved reference metadata to: {}", metadata_path.display());
        
        // Create DRAGMAP-style configuration file
        let config_path = output_dir.join("hash_table.cfg");
        let mut config_file = std::fs::File::create(&config_path)?;
        
        // Write basic configuration
        writeln!(config_file, "# Hash table configuration generated by NARFMAP")?;
        writeln!(config_file, "# Compatible with DRAGMAP hash table format")?;
        writeln!(config_file, "")?;
        writeln!(config_file, "reference_source     = '{}'", reference_path.display())?;
        writeln!(config_file, "hash_table           = '{}'", hash_table_path.display())?;
        writeln!(config_file, "pri_seed_bases       = {}", self.config.seed_len)?;
        writeln!(config_file, "max_seed_bases       = {}", self.config.max_ext_seed_len)?;
        writeln!(config_file, "ref_seed_interval    = {}", self.config.ref_seed_interval)?;
        writeln!(config_file, "max_seed_freq        = {}", self.config.max_seed_freq)?;
        writeln!(config_file, "target_seed_freq     = {}", self.config.target_seed_freq)?;
        writeln!(config_file, "pri_crc_bits         = {}", self.config.crc_primary)?;
        writeln!(config_file, "num_threads          = {}", self.config.num_threads)?;
        writeln!(config_file, "")?;
        
        // Write reference sequence information
        ref_metadata.write_dragmap_config(&mut config_file)?;
        
        info!("Saved hash table configuration to: {}", config_path.display());
        
        // Create configuration for the reference hashtable
        let ref_config = RefHashtableConfig::new(
            self.config.seed_len,
            HashTableType::Normal,
            reference_path.to_path_buf(),
            output_dir.to_path_buf(),
            self.config.num_threads,
        );

        // Create hash table data (using empty in-memory for now)
        let data = HashtableData::InMemory(Vec::new());

        // Create and return the hash table
        Ok(RefHashtable::new(ref_config, data, None))
    }
}

/// Hash table query interface for alignment
pub struct HashTableQuery {
    _hashtable: RefHashtable,
    config: HashTableConfig,
    kmer_hasher: KmerHasher,
    buckets: Vec<Bucket>,
    num_buckets: usize,
    pub reference_metadata: Option<ReferenceMetadata>,
}

impl HashTableQuery {
    /// Create a new hash table query interface with loaded buckets
    pub fn new(hashtable: RefHashtable, config: HashTableConfig, buckets: Vec<Bucket>) -> Result<Self> {
        let kmer_hasher = KmerHasher::with_dragmap_defaults(config.seed_len)?;
        let num_buckets = buckets.len();
        Ok(Self { 
            _hashtable: hashtable, 
            config, 
            kmer_hasher, 
            buckets, 
            num_buckets,
            reference_metadata: None,
        })
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
        
        // Load reference metadata if available
        let metadata_path = dir_path.join("reference_metadata.json");
        let reference_metadata = if metadata_path.exists() {
            match ReferenceMetadata::load(&metadata_path) {
                Ok(metadata) => {
                    info!("Loaded reference metadata: {} sequences", metadata.sequences.len());
                    Some(metadata)
                }
                Err(e) => {
                    warn!("Failed to load reference metadata: {}", e);
                    None
                }
            }
        } else {
            None
        };
        
        // Create hash table with loaded data
        let data = HashtableData::InMemory(Vec::new()); // We use buckets instead
        let hashtable = RefHashtable::new(ref_config, data, None);

        let mut query = Self::new(hashtable, config, buckets)?;
        query.reference_metadata = reference_metadata;
        Ok(query)
    }

    /// Load buckets from the binary hash table file
    fn load_buckets_from_file<P: AsRef<Path>>(bin_path: P, _config: &HashTableConfig) -> Result<Vec<Bucket>> {
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

    /// Query the hash table for potential positions of a k-mer with chain following
    pub fn query_kmer(&self, kmer: &[u8]) -> Result<Vec<u64>> {
        let kmer_str = String::from_utf8_lossy(kmer);
        debug!("=== HASH TABLE QUERY ===");
        debug!("Query k-mer: {} (length: {})", kmer_str, kmer.len());
        debug!("Configured k-mer size: {}", self.config.seed_len);
        
        if kmer.len() != self.config.seed_len {
            warn!("K-mer length {} does not match expected length {}", kmer.len(), self.config.seed_len);
            return Err(anyhow!("K-mer length {} does not match expected length {}", 
                              kmer.len(), self.config.seed_len));
        }

        // Convert k-mer to 2-bit representation and hash with CRC
        let kmer_2bit = self.kmer_hasher.sequence_to_2bit(kmer)?;
        debug!("K-mer 2-bit encoding: {:016x}", kmer_2bit);
        
        let crc_hash = self.kmer_hasher.hash_kmer(kmer_2bit)?;
        debug!("CRC hash: {:016x}", crc_hash);

        // Calculate bucket index using the same logic as in build_bucket_table
        let squeeze_factor = 1u64;
        let virtual_addr = self.kmer_hasher.get_address_from_hash(crc_hash, squeeze_factor);
        debug!("Virtual address: {:016x}", virtual_addr);
        
        let bucket_idx = (virtual_addr >> HashtableTraits::HASH_BUCKET_BYTES_LOG2) as usize % self.num_buckets;
        debug!("Bucket index: {} (of {} total buckets)", bucket_idx, self.num_buckets);
        
        // Extract the hash bits we're looking for (bits 35-57, same as in builder)
        let target_hash_bits = ((crc_hash >> 35) & 0x7FFFFF) as u32;
        debug!("Target hash bits: {:06x} (from CRC bits 35-57)", target_hash_bits);
        
        // Query with sophisticated DRAGMAP-compatible collision resolution
        self.query_kmer_with_probing_and_chaining(bucket_idx, target_hash_bits)
    }
    
    /// Enhanced k-mer query with block-constrained probing and chain following
    fn query_kmer_with_probing_and_chaining(&self, start_bucket_idx: usize, target_hash_bits: u32) -> Result<Vec<u64>> {
        let mut positions = Vec::new();
        let mut visited_buckets = std::collections::HashSet::new();
        
        // Calculate block constraints for probing
        let buckets_per_block = self.get_buckets_per_block();
        let block_start = (start_bucket_idx / buckets_per_block) * buckets_per_block;
        let block_end = (block_start + buckets_per_block).min(self.num_buckets);
        
        debug!("Querying with block constraints [{}, {}), start_idx: {}", block_start, block_end, start_bucket_idx);
        
        // First, try block-constrained probing (same order as insertion)
        let max_probes = 8.min(block_end - block_start);
        
        for probe in 0..max_probes {
            let bucket_idx = block_start + ((start_bucket_idx - block_start + probe) % (block_end - block_start));
            
            if visited_buckets.contains(&bucket_idx) {
                continue;
            }
            visited_buckets.insert(bucket_idx);
            
            if bucket_idx < self.buckets.len() {
                let bucket = &self.buckets[bucket_idx];
                debug!("Probing bucket {} (probe {})", bucket_idx, probe);
                
                // Search bucket and collect chain pointers
                let (bucket_positions, chain_pointers) = self.search_bucket_with_chains(bucket, target_hash_bits);
                positions.extend(bucket_positions);
                
                // Follow any chains found in this bucket
                for chain_pointer in chain_pointers {
                    let chain_bucket_idx = (bucket_idx + chain_pointer as usize) % self.num_buckets;
                    debug!("Following chain from bucket {} to bucket {}", bucket_idx, chain_bucket_idx);
                    
                    if !visited_buckets.contains(&chain_bucket_idx) {
                        visited_buckets.insert(chain_bucket_idx);
                        
                        if let Some(chain_bucket) = self.buckets.get(chain_bucket_idx) {
                            let (chain_positions, _) = self.search_bucket_with_chains(chain_bucket, target_hash_bits);
                            positions.extend(chain_positions);
                        }
                    }
                }
            }
        }
        
        debug!("Query completed: found {} positions in {} buckets", positions.len(), visited_buckets.len());
        Ok(positions)
    }
    
    /// Search a single bucket for matching records and collect chain pointers
    fn search_bucket_with_chains(&self, bucket: &Bucket, target_hash_bits: u32) -> (Vec<u64>, Vec<u32>) {
        let mut positions = Vec::new();
        let mut chain_pointers = Vec::new();
        
        for i in 0..8 {
            if let Some(record) = bucket.get(i) {
                if record.is_empty() {
                    continue;
                }
                
                debug!("  Record {}: type={:?}, hash_bits={:06x}", 
                       i, record.record_type(), record.hash_bits());
                
                match record.record_type() {
                    crate::hashtable::hash_record::RecordType::Hit => {
                        // Check for hash match
                        if record.hash_bits() == target_hash_bits {
                            if let Some(pos) = record.reference_position() {
                                debug!("    ✓ Hash match! Position: {}", pos);
                                positions.push(pos as u64);
                            }
                        }
                    }
                    crate::hashtable::hash_record::RecordType::HiFreq => {
                        // Check for hash match in high-frequency record
                        if record.hash_bits() == target_hash_bits {
                            if let Some(freq) = record.frequency() {
                                debug!("    ✓ High-frequency hash match! Frequency: {}", freq);
                                // TODO: Handle high-frequency k-mer lookups in extend table
                                // For now, we skip high-frequency k-mers
                            }
                        }
                    }
                    crate::hashtable::hash_record::RecordType::ChainBegMask |
                    crate::hashtable::hash_record::RecordType::ChainBegList => {
                        // Found a chain beginning - follow it
                        if let Some(chain_pointer) = record.chain_pointer() {
                            debug!("    → Chain beginning: pointer = {}", chain_pointer);
                            chain_pointers.push(chain_pointer);
                        }
                    }
                    crate::hashtable::hash_record::RecordType::ChainConMask |
                    crate::hashtable::hash_record::RecordType::ChainConList => {
                        // Found a chain continuation
                        if let Some(chain_pointer) = record.chain_pointer() {
                            debug!("    → Chain continuation: pointer = {}", chain_pointer);
                            chain_pointers.push(chain_pointer);
                        }
                    }
                    _ => {
                        debug!("    Unhandled record type: {:?}", record.record_type());
                    }
                }
            }
        }
        
        (positions, chain_pointers)
    }
    
    /// Calculate buckets per block for query (matches insertion logic)
    fn get_buckets_per_block(&self) -> usize {
        // Use the same block size calculation as insertion for consistency
        let block_size_log2 = if self.config.crc_primary >= 16 {
            // For larger hash tables, use 256-bucket blocks (2^8 buckets = 2^14 bytes)
            8
        } else if self.config.crc_primary >= 12 {
            // For medium hash tables, use 64-bucket blocks (2^6 buckets = 2^12 bytes)
            6
        } else {
            // For small hash tables, use 16-bucket blocks (2^4 buckets = 2^10 bytes)
            4
        };
        
        let buckets_per_block_log2 = block_size_log2.max(3); // Minimum 8 buckets per block
        1 << buckets_per_block_log2
    }

    /// Get hash table statistics
    pub fn get_stats(&self) -> HashMap<String, u64> {
        let mut stats = HashMap::new();
        // stats.insert("size_bytes".to_string(), self._hashtable.size_bytes() as u64);
        stats.insert("size_bytes".to_string(), 0u64); // TODO: Calculate actual size
        stats.insert("kmer_size".to_string(), self.config.seed_len as u64);
        stats.insert("max_seed_freq".to_string(), self.config.max_seed_freq as u64);
        stats
    }
}

/// Legacy Hashtable struct for compatibility
/// This is kept for backward compatibility with existing code
#[deprecated(note = "Use HashTableBuilder and HashTableQuery instead")]
pub struct Hashtable {
    _query: Option<HashTableQuery>,
}

#[allow(deprecated)]
impl Hashtable {
    pub fn new() -> Self {
        Self { _query: None }
    }

    /// Create from a config (for compatibility)
    pub fn from_config(_config: HashTableConfig) -> Self {
        Self { _query: None }
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
        
        let mut builder = HashTableBuilder::new(config).unwrap();
        
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
        
        let mut builder = HashTableBuilder::new(config).unwrap();
        
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
            
            let mut builder = HashTableBuilder::new(config).unwrap();
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
        
        let mut builder = HashTableBuilder::new(config).unwrap();
        
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
    
    #[test]
    fn test_linear_probing() {
        let mut config = HashTableConfig::default();
        config.seed_len = 16;
        config.crc_primary = 10; // Small table for testing (2^10 bits = 16 buckets)
        
        let mut builder = HashTableBuilder::new(config).unwrap();
        
        // Create buckets
        let mut buckets = vec![Bucket::new(); 16];
        
        // Fill bucket 0 completely
        for i in 0..8 {
            let record = HashRecord::hit(0, i as u32, false, false, false, i as u32 * 100);
            buckets[0].set(i, record).unwrap();
        }
        assert!(buckets[0].is_full());
        
        // Try to insert another record that would go to bucket 0
        let new_record = HashRecord::hit(1, 0x9999, false, false, false, 9999);
        let inserted = builder.insert_record_with_probing(&mut buckets, 0, new_record);
        
        assert!(inserted);
        // Should have been placed in bucket 1
        assert_eq!(buckets[1].count_occupied(), 1);
        assert_eq!(buckets[1][0].reference_position(), Some(9999));
    }
    
    #[test]
    fn test_reference_metadata_preservation() {
        use tempfile::TempDir;
        
        // Create test sequences with meaningful names
        let sequences = vec![
            crate::io::sequence::Sequence::new(
                "chr1_human_genome".to_string(),
                "ACGTACGTACGTACGTACGTACGTACGTACGT",
                None
            ).unwrap(),
            crate::io::sequence::Sequence::new(
                "chr2_human_genome".to_string(),
                "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA",
                None
            ).unwrap(),
            crate::io::sequence::Sequence::new(
                "chrM_mitochondrial".to_string(),
                "TTTTAAAACCCCGGGGTTTTAAAACCCCGGGG",
                None
            ).unwrap(),
        ];
        
        // Create metadata from sequences
        let metadata = ReferenceMetadata::from_sequences(&sequences);
        
        // Verify sequence names are preserved
        assert_eq!(metadata.sequences.len(), 3);
        assert_eq!(metadata.sequences[0].name, "chr1_human_genome");
        assert_eq!(metadata.sequences[1].name, "chr2_human_genome");
        assert_eq!(metadata.sequences[2].name, "chrM_mitochondrial");
        
        // Verify positions and lengths
        assert_eq!(metadata.sequences[0].index, 0);
        assert_eq!(metadata.sequences[0].length, 32);
        assert_eq!(metadata.sequences[0].start_position, 0);
        
        assert_eq!(metadata.sequences[1].index, 1);
        assert_eq!(metadata.sequences[1].length, 32);
        // Should have padding after first sequence
        assert!(metadata.sequences[1].start_position > 32);
        
        // Test save/load roundtrip
        let temp_dir = TempDir::new().unwrap();
        let metadata_path = temp_dir.path().join("test_metadata.json");
        
        metadata.save(&metadata_path).unwrap();
        let loaded_metadata = ReferenceMetadata::load(&metadata_path).unwrap();
        
        assert_eq!(loaded_metadata.sequences.len(), metadata.sequences.len());
        for i in 0..metadata.sequences.len() {
            assert_eq!(loaded_metadata.sequences[i].name, metadata.sequences[i].name);
            assert_eq!(loaded_metadata.sequences[i].length, metadata.sequences[i].length);
        }
    }
    
    #[test]
    fn test_collision_handling_and_chaining() {
        let mut config = HashTableConfig::default();
        config.seed_len = 16;
        config.crc_primary = 12; // Small table to force collisions (2^12 bits = 64 buckets)
        config.max_seed_freq = 100;
        
        let mut builder = HashTableBuilder::new(config.clone()).unwrap();
        
        // Create many k-mer entries that will force collisions
        let mut entries = HashMap::new();
        
        // Add 200 different k-mers to a table with only 64 buckets
        // This will definitely cause collisions and test our probing/chaining
        for i in 0..200 {
            let hash_value = (i as u64) << 20; // Spread them across hash space
            entries.insert(hash_value, vec![i as u64 * 100]); // Position = i * 100
        }
        
        // Build bucket table with collision handling
        let buckets = builder.build_bucket_table(entries).unwrap();
        
        // Should have exactly 64 buckets (2^(12-6))
        assert_eq!(buckets.len(), 64);
        
        // Count total records inserted
        let total_records: usize = buckets.iter()
            .map(|b| b.count_occupied())
            .sum();
        
        // Should have inserted 200 records
        assert_eq!(total_records, 200);
        
        // Check that some buckets have multiple records (due to probing)
        let multi_record_buckets: usize = buckets.iter()
            .filter(|b| b.count_occupied() > 1)
            .count();
        assert!(multi_record_buckets > 0, "Expected some buckets to have multiple records due to collisions");
        
        // Print collision statistics
        let stats = &builder.collision_stats;
        println!("=== Collision Test Results ===");
        println!("Total insertions: {}", stats.total_insertions);
        println!("Direct insertions: {} ({:.1}%)", 
                 stats.direct_insertions, 
                 stats.direct_insertions as f64 / stats.total_insertions as f64 * 100.0);
        println!("Probing insertions: {} ({:.1}%)", 
                 stats.probing_insertions,
                 stats.probing_insertions as f64 / stats.total_insertions as f64 * 100.0);
        println!("Chaining insertions: {} ({:.1}%)", 
                 stats.chaining_insertions,
                 stats.chaining_insertions as f64 / stats.total_insertions as f64 * 100.0);
        println!("Failed insertions: {}", stats.failed_insertions);
        println!("Max probe depth: {}", stats.max_probe_depth);
        println!("Total chains: {}", stats.total_chains);
        
        // Verify collision handling worked
        assert!(stats.total_insertions == 200);
        assert!(stats.probing_insertions > 0 || stats.chaining_insertions > 0, 
                "Expected some collisions to be resolved by probing or chaining");
        assert_eq!(stats.failed_insertions, 0, "No insertions should have failed");
        
        // Test querying with collision resolution
        let ref_config = crate::reference::hashtable::HashtableConfig::new(
            config.seed_len,
            crate::reference::hashtable::HashTableType::Normal,
            std::path::PathBuf::from("test"),
            std::path::PathBuf::from("test_out"),
            config.num_threads,
        );
        
        let query = HashTableQuery::new(
            crate::reference::hashtable::Hashtable::new(
                ref_config,
                crate::reference::hashtable::HashtableData::InMemory(Vec::new()),
                None
            ),
            config,
            buckets
        ).unwrap();
        
        // Test that we can find records even with collisions
        // Use a simple k-mer for testing
        let test_kmer = b"ACGTACGTACGTACGT"; // 16 bases
        
        // This test mainly verifies that the query doesn't crash with the new collision logic
        let result = query.query_kmer(test_kmer);
        assert!(result.is_ok(), "Query should succeed even if no matches found");
    }
    
    #[test]
    fn test_hash_table_with_metadata_integration() {
        use tempfile::TempDir;
        use std::fs;
        
        // Skip if no test data available
        let test_fasta = Path::new("data/tiny/tiny-2x1Xrepeats.v8/tiny.fasta");
        if !test_fasta.exists() {
            println!("Test FASTA not found, skipping integration test");
            return;
        }
        
        let temp_dir = TempDir::new().unwrap();
        let output_dir = temp_dir.path();
        
        // Configure hash table builder
        let mut config = HashTableConfig::default();
        config.seed_len = 21;
        config.ref_seed_interval = 10.0;
        config.crc_primary = 16;
        config.num_threads = 1;
        
        // Build hash table
        let mut builder = HashTableBuilder::new(config.clone()).unwrap();
        builder.build_from_fasta(test_fasta, output_dir).unwrap();
        
        // Verify metadata file was created
        let metadata_path = output_dir.join("reference_metadata.json");
        assert!(metadata_path.exists());
        
        // Load and verify metadata
        let metadata = ReferenceMetadata::load(&metadata_path).unwrap();
        assert!(metadata.sequences.len() > 0);
        assert!(metadata.sequences[0].name.contains("phiX174"));
        
        // Verify DRAGMAP config file was created
        let config_path = output_dir.join("hash_table.cfg");
        assert!(config_path.exists());
        
        // Read config file and verify it contains sequence information
        let config_content = fs::read_to_string(&config_path).unwrap();
        assert!(config_content.contains("reference_sequence0"));
        assert!(config_content.contains("phiX174"));
        
        // Load hash table and verify metadata is loaded
        let query = HashTableQuery::load_from_dir(output_dir, config).unwrap();
        assert!(query.reference_metadata.is_some());
        
        let loaded_metadata = query.reference_metadata.as_ref().unwrap();
        assert_eq!(loaded_metadata.sequences.len(), metadata.sequences.len());
        assert_eq!(loaded_metadata.sequences[0].name, metadata.sequences[0].name);
    }
} 