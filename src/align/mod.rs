use std::path::Path;
use std::collections::HashMap;
use anyhow::{Result, anyhow};
use log::{info, debug, warn};

use crate::config::{Config, AlignmentConfig, HashTableConfig};
use crate::hashtable::{HashTableQuery, KmerHasher};
use crate::io::fastq::FastqReader;
use crate::io::sequence::Sequence;
use crate::reference::sequence::NucSeq;

pub mod sam;

/// Alignment result for a single read
#[derive(Debug, Clone)]
pub struct AlignmentResult {
    pub read_id: String,
    pub reference_id: String,
    pub position: u32,
    pub cigar: String,
    pub mapq: u8,
    pub is_reverse: bool,
    pub is_paired: bool,
    pub is_proper_pair: bool,
    pub mate_reference_id: Option<String>,
    pub mate_position: Option<u32>,
    pub template_length: Option<i32>,
}

/// Seed hit from hash table lookup
#[derive(Debug, Clone)]
pub struct SeedHit {
    pub read_position: usize,
    pub reference_position: u32,
    pub sequence_id: u64,
    pub is_reverse: bool,
}

/// Main aligner implementing seed-and-extend algorithm
pub struct Aligner {
    config: AlignmentConfig,
    hash_table: HashTableQuery,
    kmer_hasher: KmerHasher,
}

impl Aligner {
    /// Create a new aligner with hash table and configuration
    pub fn new(config: AlignmentConfig, hash_table_dir: &Path) -> Result<Self> {
        info!("Loading hash table from: {}", hash_table_dir.display());
        
        let hash_table_config = HashTableConfig::default(); // TODO: Load from config file
        let hash_table = HashTableQuery::load_from_dir(hash_table_dir, hash_table_config.clone())?;
        
        let kmer_hasher = KmerHasher::with_dragmap_defaults(hash_table_config.seed_len)?;
        
        Ok(Self {
            config,
            hash_table,
            kmer_hasher,
        })
    }

    /// Align a single read to the reference
    pub fn align_read(&self, read: &Sequence) -> Result<Option<AlignmentResult>> {
        debug!("Aligning read: {} ({} bp)", read.id, read.len());
        
        // Extract seeds from the read
        let seeds = self.extract_seeds(read)?;
        debug!("Extracted {} seeds from read", seeds.len());
        
        if seeds.is_empty() {
            warn!("No valid seeds found in read {}", read.id);
            return Ok(None);
        }
        
        // Find seed hits in the hash table
        let mut all_hits = Vec::new();
        for (position, seed_seq) in seeds {
            if let Ok(hits) = self.find_seed_hits(&seed_seq, position) {
                all_hits.extend(hits);
            }
        }
        
        debug!("Found {} seed hits for read {}", all_hits.len(), read.id);
        
        if all_hits.is_empty() {
            debug!("No seed hits found for read {}", read.id);
            return Ok(None);
        }
        
        // Cluster hits by reference position
        let clustered_hits = self.cluster_hits(all_hits);
        debug!("Clustered into {} potential alignment positions", clustered_hits.len());
        
        // Score and select best alignment
        if let Some(best_hit) = self.select_best_alignment(read, clustered_hits)? {
            debug!("Best alignment for {}: ref_pos={}, score={}", read.id, best_hit.reference_position, "TODO");
            
            // TODO: Implement full alignment extension and CIGAR generation
            let alignment = AlignmentResult {
                read_id: read.id.clone(),
                reference_id: format!("ref_{}", best_hit.sequence_id), // TODO: Get actual reference name
                position: best_hit.reference_position,
                cigar: format!("{}M", read.len()), // TODO: Generate proper CIGAR
                mapq: 60, // TODO: Calculate mapping quality
                is_reverse: best_hit.is_reverse,
                is_paired: false, // TODO: Handle paired-end
                is_proper_pair: false,
                mate_reference_id: None,
                mate_position: None,
                template_length: None,
            };
            
            Ok(Some(alignment))
        } else {
            debug!("No valid alignment found for read {}", read.id);
            Ok(None)
        }
    }
    
    /// Extract k-mer seeds from a read
    fn extract_seeds(&self, read: &Sequence) -> Result<Vec<(usize, Vec<u8>)>> {
        let mut seeds = Vec::new();
        let k = self.config.seed_len;
        
        if read.len() < k {
            return Ok(seeds);
        }
        
        // Extract seeds with configurable step size
        let step = self.config.seed_step_size.max(1);
        
        for i in (0..=read.len() - k).step_by(step) {
            if let Some(kmer_seq) = read.substring(i, i + k) {
                // Convert to raw sequence bytes
                let mut kmer_bytes = Vec::with_capacity(k);
                for j in 0..k {
                    kmer_bytes.push(match kmer_seq.get(j) {
                        crate::reference::sequence::Nucleotide::A => b'A',
                        crate::reference::sequence::Nucleotide::C => b'C',
                        crate::reference::sequence::Nucleotide::G => b'G',
                        crate::reference::sequence::Nucleotide::T => b'T',
                        crate::reference::sequence::Nucleotide::N => b'N',
                    });
                }
                seeds.push((i, kmer_bytes));
            }
        }
        
        Ok(seeds)
    }
    
    /// Find hits for a k-mer seed in the hash table
    fn find_seed_hits(&self, seed: &[u8], read_position: usize) -> Result<Vec<SeedHit>> {
        let positions = self.hash_table.query_kmer(seed)?;
        
        let mut hits = Vec::new();
        for encoded_pos in positions {
            let sequence_id = encoded_pos >> 32;
            let reference_position = (encoded_pos & 0xFFFFFFFF) as u32;
            
            hits.push(SeedHit {
                read_position,
                reference_position,
                sequence_id,
                is_reverse: false, // TODO: Detect reverse complement
            });
        }
        
        Ok(hits)
    }
    
    /// Cluster hits that are close together on the reference
    fn cluster_hits(&self, hits: Vec<SeedHit>) -> Vec<Vec<SeedHit>> {
        let mut clusters: Vec<Vec<SeedHit>> = Vec::new();
        let cluster_distance = self.config.cluster_distance;
        
        for hit in hits {
            let mut found_cluster = false;
            
            for cluster in clusters.iter_mut() {
                if let Some(first_hit) = cluster.first() {
                    // Check if this hit is close to the cluster
                    if hit.sequence_id == first_hit.sequence_id &&
                       hit.reference_position.abs_diff(first_hit.reference_position) <= cluster_distance {
                        cluster.push(hit.clone());
                        found_cluster = true;
                        break;
                    }
                }
            }
            
            if !found_cluster {
                clusters.push(vec![hit]);
            }
        }
        
        // Sort clusters by number of hits (descending)
        clusters.sort_by(|a, b| b.len().cmp(&a.len()));
        
        clusters
    }
    
    /// Select the best alignment from clustered hits
    fn select_best_alignment(&self, read: &Sequence, mut clusters: Vec<Vec<SeedHit>>) -> Result<Option<SeedHit>> {
        if clusters.is_empty() {
            return Ok(None);
        }
        
        // For now, just take the cluster with the most hits
        let best_cluster = clusters.remove(0);
        
        // Within the cluster, take the hit with the most consistent position
        let best_hit = best_cluster.into_iter()
            .min_by_key(|hit| hit.reference_position)
            .unwrap();
        
        Ok(Some(best_hit))
    }
    
    /// Align paired-end reads
    pub fn align_paired_reads(&self, read1: &Sequence, read2: &Sequence) -> Result<(Option<AlignmentResult>, Option<AlignmentResult>)> {
        debug!("Aligning paired reads: {} and {}", read1.id, read2.id);
        
        // For now, align each read independently
        let alignment1 = self.align_read(read1)?;
        let alignment2 = self.align_read(read2)?;
        
        // TODO: Implement proper paired-end alignment logic
        // - Ensure reads are properly paired
        // - Calculate insert size
        // - Set proper SAM flags
        
        Ok((alignment1, alignment2))
    }
} 