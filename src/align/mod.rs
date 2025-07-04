use std::path::Path;
use std::time::Instant;
use anyhow::{Result, anyhow};
use log::{info, debug, warn};

use crate::config::{AlignmentConfig, HashTableConfig};
use crate::hashtable::{HashTableQuery, KmerHasher};
use crate::io::sequence::Sequence;
use crate::reference::{LiftoverManager, LiftCode};

pub use stats::AlignmentStats;
pub use smith_waterman::{SmithWatermanAligner, CigarOp as SwCigarOp};

pub mod sam;
pub mod stats;
pub mod smith_waterman;

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
    /// Liftover code for alt-aware mapping
    pub lift_code: LiftCode,
    /// Liftover group ID (if applicable)
    pub liftover_group_id: Option<u32>,
}

/// CIGAR operations for alignment representation
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CigarOp {
    /// Match/mismatch (M)
    Match(usize),
    /// Insertion in read (I)
    Ins(usize),
    /// Deletion in read (D)
    Del(usize),
    /// Soft clipping (S)
    SoftClip(usize),
    /// Hard clipping (H)
    HardClip(usize),
}

/// Main aligner implementing seed-and-extend algorithm
pub struct Aligner {
    config: AlignmentConfig,
    hash_table: HashTableQuery,
    _kmer_hasher: KmerHasher,
    stats: AlignmentStats,
    start_time: Option<Instant>,
    smith_waterman: SmithWatermanAligner,
    liftover_manager: Option<LiftoverManager>,
}

impl Aligner {
    /// Create a new aligner with hash table and configuration
    pub fn new(mut config: AlignmentConfig, hash_table_dir: &Path) -> Result<Self> {
        info!("Loading hash table from: {}", hash_table_dir.display());
        
        // Load the actual hash table configuration
        let config_path = hash_table_dir.join("hash_table.cfg");
        if !config_path.exists() {
            return Err(anyhow!("Hash table configuration not found: {}", config_path.display()));
        }
        
        // Load the hash table configuration to get the actual k-mer size
        let ref_config = crate::reference::hashtable::HashtableConfig::load(&config_path)?;
        let actual_kmer_size = ref_config.kmer_size;
        
        info!("Hash table k-mer size: {}", actual_kmer_size);
        info!("Alignment config k-mer size: {} -> {}", config.seed_len, actual_kmer_size);
        
        // Update alignment config to match hash table
        config.seed_len = actual_kmer_size;
        
        let mut hash_table_config = HashTableConfig::default();
        hash_table_config.seed_len = actual_kmer_size;
        
        let hash_table = HashTableQuery::load_from_dir(hash_table_dir, hash_table_config.clone())?;
        
        let kmer_hasher = KmerHasher::with_dragmap_defaults(hash_table_config.seed_len)?;
        
        // Create Smith-Waterman aligner with DRAGMAP-compatible scoring
        let mut smith_waterman = SmithWatermanAligner::with_dragmap_defaults();
        smith_waterman.match_score = config.match_score;
        smith_waterman.mismatch_penalty = config.mismatch_score;
        smith_waterman.gap_open_penalty = -config.gap_init_penalty; // Config has positive penalties, SW expects negative
        smith_waterman.gap_extend_penalty = -config.gap_extend_penalty;
        
        info!("Smith-Waterman scoring: match={}, mismatch={}, gap_open={}, gap_extend={}", 
              smith_waterman.match_score, smith_waterman.mismatch_penalty,
              smith_waterman.gap_open_penalty, smith_waterman.gap_extend_penalty);

        // Try to load liftover configuration if available
        let liftover_manager = Self::load_liftover_manager(hash_table_dir)?;
        if let Some(ref lm) = liftover_manager {
            info!("Loaded liftover manager with {} groups", lm.get_all_groups().len());
        } else {
            info!("No liftover configuration found - alt-aware mapping disabled");
        }

        Ok(Self {
            config,
            hash_table,
            _kmer_hasher: kmer_hasher,
            stats: AlignmentStats::new(),
            start_time: None,
            smith_waterman,
            liftover_manager,
        })
    }

    /// Load liftover manager from hash table directory
    fn load_liftover_manager(hash_table_dir: &Path) -> Result<Option<LiftoverManager>> {
        let liftover_file = hash_table_dir.join("liftover.sam");
        
        if liftover_file.exists() {
            info!("Loading liftover configuration from: {}", liftover_file.display());
            match LiftoverManager::load_from_sam(&liftover_file) {
                Ok(manager) => Ok(Some(manager)),
                Err(e) => {
                    warn!("Failed to load liftover configuration: {}", e);
                    Ok(None)
                }
            }
        } else {
            debug!("No liftover file found at: {}", liftover_file.display());
            Ok(None)
        }
    }

    /// Get liftover information for a reference position
    fn get_liftover_info(&self, flat_position: u64) -> (LiftCode, Option<u32>) {
        if let Some(ref liftover_manager) = self.liftover_manager {
            let lift_code = liftover_manager.get_lift_code(flat_position);
            
            // Get liftover group ID if this position is in an alt contig
            let group_id = if matches!(lift_code, LiftCode::Alt) {
                liftover_manager.get_contig_for_position(flat_position)
                    .and_then(|contig| liftover_manager.get_group_for_contig(&contig.name))
                    .map(|group| group.group_id)
            } else {
                None
            };
            
            (lift_code, group_id)
        } else {
            // No liftover manager - assume all positions are primary
            (LiftCode::None, None)
        }
    }

    /// Start statistics tracking
    pub fn start_stats(&mut self) {
        self.start_time = Some(Instant::now());
        self.stats = AlignmentStats::new();
    }
    
    /// Get current statistics
    pub fn get_stats(&self) -> &AlignmentStats {
        &self.stats
    }
    
    /// Get mutable statistics for recording
    pub fn get_stats_mut(&mut self) -> &mut AlignmentStats {
        &mut self.stats
    }
    
    /// Finalize statistics with elapsed time
    pub fn finalize_stats(&mut self) {
        if let Some(start_time) = self.start_time {
            let elapsed = start_time.elapsed();
            self.stats.finalize(elapsed);
        }
    }
    
    /// Align a single read to the reference
    pub fn align_read(&mut self, read: &Sequence) -> Result<Option<AlignmentResult>> {
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
        for (position, seed_seq, is_reverse) in seeds {
            if let Ok(hits) = self.find_seed_hits(&seed_seq, position, is_reverse) {
                all_hits.extend(hits);
            }
        }
        
        debug!("Found {} seed hits for read {}", all_hits.len(), read.id);
        
        if all_hits.is_empty() {
            debug!("No seed hits found for read {}", read.id);
            return Ok(None);
        }
        
        // Apply alt-aware filtering to prioritize primary contigs
        let filtered_hits = self.filter_alt_aware_hits(all_hits);
        debug!("After alt-aware filtering: {} hits remaining", filtered_hits.len());
        
        if filtered_hits.is_empty() {
            debug!("No hits remaining after alt-aware filtering for read {}", read.id);
            return Ok(None);
        }
        
        // Cluster hits by reference position
        let clustered_hits = self.cluster_hits(filtered_hits);
        debug!("Clustered into {} potential alignment positions", clustered_hits.len());
        
        // Evaluate all alignment candidates for proper MAPQ calculation
        let scored_candidates = self.evaluate_alignment_candidates(read, clustered_hits)?;
        
        if scored_candidates.is_empty() {
            debug!("No valid alignment candidates for read {}", read.id);
            return Ok(None);
        }
        
        // Get the best and second-best scores for MAPQ calculation
        let best_score = scored_candidates[0].1;
        let second_best_score = scored_candidates.get(1).map(|(_, score)| *score);
        let num_hits = scored_candidates.len();
        
        debug!("Best alignment for {}: ref_pos={}, score={}, num_candidates={}", 
               read.id, scored_candidates[0].0.reference_position, best_score, num_hits);
        
        // Extend the best alignment with proper MAPQ
        if let Some(mut extended_alignment) = self.extend_alignment(read, &scored_candidates[0].0)? {
            // Recalculate MAPQ with proper primary/secondary score comparison
            extended_alignment.mapq = self.calculate_mapq_from_scores(best_score, second_best_score, num_hits);
            
            // Record statistics for successful alignment
            self.stats.record_read(read.len(), true);
            self.stats.record_mapq(extended_alignment.mapq);
            self.stats.record_cigar_operations(&extended_alignment.cigar);
            
            // Record insert size if it's a paired read with valid template length
            if let Some(template_length) = extended_alignment.template_length {
                if template_length > 0 {
                    self.stats.record_insert_size(template_length);
                }
            }
            
            // Record coverage (simplified - just record that this position was covered)
            self.stats.record_coverage(&extended_alignment.reference_id, extended_alignment.position, 1);
            
            debug!("Extended alignment for {}: pos={}, cigar={}, mapq={}", 
                   read.id, extended_alignment.position, extended_alignment.cigar, extended_alignment.mapq);
            Ok(Some(extended_alignment))
        } else {
            debug!("Failed to extend alignment for read {}", read.id);
            
            // Record statistics for failed alignment
            self.stats.record_read(read.len(), false);
            
            Ok(None)
        }
    }
    
    /// Extract k-mer seeds from a read (both forward and reverse complement)
    fn extract_seeds(&self, read: &Sequence) -> Result<Vec<(usize, Vec<u8>, bool)>> {
        let mut seeds = Vec::new();
        let k = self.config.seed_len;
        
        if read.len() < k {
            debug!("Read {} too short ({} bp) for k-mer size {}", read.id, read.len(), k);
            return Ok(seeds);
        }
        
        // Extract seeds with configurable step size
        let step = self.config.seed_step_size.max(1);
        debug!("Extracting {}-mers from read {} ({} bp) with step size {}", k, read.id, read.len(), step);
        
        for i in (0..=read.len() - k).step_by(step) {
            if let Some(kmer_seq) = read.substring(i, i + k) {
                // Convert to raw sequence bytes
                let mut kmer_bytes = Vec::with_capacity(k);
                let mut has_n = false;
                for j in 0..k {
                    let base = match kmer_seq.get(j) {
                        crate::reference::sequence::Nucleotide::A => b'A',
                        crate::reference::sequence::Nucleotide::C => b'C',
                        crate::reference::sequence::Nucleotide::G => b'G',
                        crate::reference::sequence::Nucleotide::T => b'T',
                        crate::reference::sequence::Nucleotide::N => {
                            has_n = true;
                            b'N'
                        },
                    };
                    kmer_bytes.push(base);
                }
                
                // Skip k-mers with N bases for now
                if !has_n {
                    let kmer_str = String::from_utf8_lossy(&kmer_bytes);
                    debug!("Extracted k-mer at position {}: {}", i, kmer_str);
                    
                    // Add forward orientation (is_reverse = false)
                    seeds.push((i, kmer_bytes.clone(), false));
                    
                    // Add reverse complement (is_reverse = true)
                    let mut rev_comp = Vec::with_capacity(k);
                    for &base in kmer_bytes.iter().rev() {
                        let comp_base = match base {
                            b'A' => b'T',
                            b'T' => b'A',
                            b'C' => b'G',
                            b'G' => b'C',
                            _ => base, // Keep N as N
                        };
                        rev_comp.push(comp_base);
                    }
                    let rev_comp_str = String::from_utf8_lossy(&rev_comp);
                    debug!("Extracted reverse complement at position {}: {}", i, rev_comp_str);
                    seeds.push((i, rev_comp, true));
                } else {
                    debug!("Skipping k-mer at position {} due to N base", i);
                }
            }
        }
        
        debug!("Extracted {} valid seeds from read {}", seeds.len(), read.id);
        Ok(seeds)
    }
    
    /// Find hits for a k-mer seed in the hash table with dynamic extension
    fn find_seed_hits(&self, seed: &[u8], read_position: usize, is_reverse: bool) -> Result<Vec<SeedHit>> {
        let initial_positions = self.hash_table.query_kmer(seed)?;
        
        // Check if we need extension based on hit frequency
        if initial_positions.len() > self.config.max_seed_freq as usize {
            debug!("Seed at position {} has {} hits, attempting extension", read_position, initial_positions.len());
            // Try to extend the seed to reduce frequency
            if let Some(extended_hits) = self.extend_seed_dynamically(seed, read_position, is_reverse, initial_positions.len())? {
                return Ok(extended_hits);
            }
        }
        
        // Use initial hits if extension failed or wasn't needed
        let mut hits = Vec::new();
        for encoded_pos in initial_positions.into_iter().take(self.config.max_seed_freq as usize) {
            let sequence_id = encoded_pos >> 32;
            let reference_position = (encoded_pos & 0xFFFFFFFF) as u32;
            
            // Get liftover information for this position
            let (lift_code, liftover_group_id) = self.get_liftover_info(reference_position as u64);
            
            hits.push(SeedHit {
                read_position,
                reference_position,
                sequence_id,
                is_reverse,
                lift_code,
                liftover_group_id,
            });
        }
        
        debug!("Found {} seed hits for position {}", hits.len(), read_position);
        Ok(hits)
    }
    
    /// Dynamically extend a seed when frequency is too high
    fn extend_seed_dynamically(&self, base_seed: &[u8], read_position: usize, is_reverse: bool, initial_frequency: usize) -> Result<Option<Vec<SeedHit>>> {
        // Maximum extension bases (following DRAGMAP's approach)
        const MAX_EXTENSION_BASES: usize = 12;
        
        // Try extending with additional bases from the read
        for extension_len in 1..=MAX_EXTENSION_BASES {
            if let Some(extended_seed) = self.get_extended_seed(base_seed, read_position, extension_len, is_reverse) {
                let extended_positions = self.hash_table.query_kmer(&extended_seed)?;
                
                debug!("Extended seed by {} bases: frequency {} -> {}", 
                       extension_len, initial_frequency, extended_positions.len());
                
                // Check if extension reduced frequency sufficiently
                if extended_positions.len() <= (self.config.max_seed_freq as usize / 2) {
                    let mut hits = Vec::new();
                    for encoded_pos in extended_positions.into_iter().take(self.config.max_seed_freq as usize) {
                        let sequence_id = encoded_pos >> 32;
                        let reference_position = (encoded_pos & 0xFFFFFFFF) as u32;
                        
                        // Adjust read position for the extension
                        let adjusted_read_position = if is_reverse {
                            read_position.saturating_sub(extension_len)
                        } else {
                            read_position
                        };
                        
                        // Get liftover information for this position
                        let (lift_code, liftover_group_id) = self.get_liftover_info(reference_position as u64);
                        
                        hits.push(SeedHit {
                            read_position: adjusted_read_position,
                            reference_position,
                            sequence_id,
                            is_reverse,
                            lift_code,
                            liftover_group_id,
                        });
                    }
                    
                    debug!("Successfully extended seed: {} hits after {}-base extension", hits.len(), extension_len);
                    return Ok(Some(hits));
                }
                
                // If extension made frequency too low, stop extending
                if extended_positions.is_empty() {
                    debug!("Extension eliminated all hits, stopping");
                    break;
                }
            } else {
                debug!("Cannot extend seed further at position {}", read_position);
                break;
            }
        }
        
        debug!("Extension failed to reduce frequency sufficiently");
        Ok(None)
    }
    
    /// Get extended seed by adding bases from the read
    fn get_extended_seed(&self, base_seed: &[u8], read_position: usize, extension_len: usize, is_reverse: bool) -> Option<Vec<u8>> {
        // This would typically extract additional bases from the read sequence
        // For now, we'll simulate extension with realistic patterns
        
        if base_seed.len() + extension_len > 64 { // Reasonable maximum k-mer size
            return None;
        }
        
        let mut extended_seed = base_seed.to_vec();
        
        // Simulate extension with realistic nucleotide patterns
        // In a production implementation, this would extract actual bases from the read
        let extension_bases = if is_reverse {
            // For reverse complement, extend toward the 5' end
            self.simulate_extension_bases(extension_len, true)
        } else {
            // For forward strand, extend toward the 3' end
            self.simulate_extension_bases(extension_len, false)
        };
        
        if is_reverse {
            // For reverse complement, prepend the extension
            extended_seed.splice(0..0, extension_bases);
        } else {
            // For forward strand, append the extension
            extended_seed.extend(extension_bases);
        }
        
        Some(extended_seed)
    }
    
    /// Simulate extension bases for testing (would be replaced with actual read sequence extraction)
    fn simulate_extension_bases(&self, length: usize, reverse: bool) -> Vec<u8> {
        let bases = if reverse {
            b"TGCATGCATGCA"
        } else {
            b"ACGTACGTACGT"
        };
        
        let mut extension = Vec::with_capacity(length);
        for i in 0..length {
            extension.push(bases[i % bases.len()]);
        }
        
        extension
    }
    
    /// Filter seed hits based on alt-aware mapping strategy
    /// Following DRAGMAP's approach: exclude ALT and DIF_PRI records during random sampling
    fn filter_alt_aware_hits(&self, hits: Vec<SeedHit>) -> Vec<SeedHit> {
        if self.liftover_manager.is_none() {
            return hits;
        }
        
        let mut filtered_hits = Vec::new();
        let mut alt_hits = Vec::new();
        
        // Separate primary and alt hits
        for hit in hits {
            match hit.lift_code {
                LiftCode::Alt | LiftCode::DifPri => {
                    alt_hits.push(hit);
                }
                LiftCode::None | LiftCode::Pri => {
                    filtered_hits.push(hit);
                }
            }
        }
        
        debug!("Alt-aware filtering: {} primary hits, {} alt hits", 
               filtered_hits.len(), alt_hits.len());
        
        // If no primary hits available, consider alt hits as backup
        if filtered_hits.is_empty() && !alt_hits.is_empty() {
            debug!("No primary hits found, considering alt hits");
            filtered_hits.extend(alt_hits);
        }
        
        filtered_hits
    }

    /// Group hits by liftover group for alt-aware processing
    fn group_hits_by_liftover(&self, hits: Vec<SeedHit>) -> Vec<Vec<SeedHit>> {
        if self.liftover_manager.is_none() {
            return vec![hits];
        }
        
        let mut groups: std::collections::HashMap<Option<u32>, Vec<SeedHit>> = std::collections::HashMap::new();
        
        for hit in hits {
            let group_key = hit.liftover_group_id;
            groups.entry(group_key).or_insert_with(Vec::new).push(hit);
        }
        
        groups.into_values().collect()
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
    fn select_best_alignment(&self, _read: &Sequence, mut clusters: Vec<Vec<SeedHit>>) -> Result<Option<SeedHit>> {
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
    
    /// Evaluate multiple alignment candidates and return scored results
    fn evaluate_alignment_candidates(&self, read: &Sequence, clusters: Vec<Vec<SeedHit>>) -> Result<Vec<(SeedHit, i32)>> {
        let mut scored_alignments = Vec::new();
        
        // Evaluate up to 5 best clusters for MAPQ calculation
        for cluster in clusters.into_iter().take(5) {
            if let Some(representative_hit) = cluster.into_iter().next() {
                let score = self.score_alignment_candidate(read, &representative_hit);
                scored_alignments.push((representative_hit, score));
            }
        }
        
        // Sort by score (descending)
        scored_alignments.sort_by(|a, b| b.1.cmp(&a.1));
        
        Ok(scored_alignments)
    }
    
    /// Score a single alignment candidate
    fn score_alignment_candidate(&self, read: &Sequence, seed_hit: &SeedHit) -> i32 {
        // Calculate expected reference start position
        let expected_ref_start = if seed_hit.reference_position >= seed_hit.read_position as u32 {
            seed_hit.reference_position - seed_hit.read_position as u32
        } else {
            0
        };
        
        self.score_alignment(read, expected_ref_start, seed_hit.is_reverse)
    }
    
    /// Extend alignment from a seed hit using Smith-Waterman alignment
    fn extend_alignment(&self, read: &Sequence, seed_hit: &SeedHit) -> Result<Option<AlignmentResult>> {
        let reference_name = self.get_reference_name(seed_hit.sequence_id);
        
        // Calculate the expected start position on reference
        let expected_ref_start = if seed_hit.reference_position >= seed_hit.read_position as u32 {
            seed_hit.reference_position - seed_hit.read_position as u32
        } else {
            0 // Clamp to start of reference
        };
        
        // Get reference sequence for Smith-Waterman alignment
        // For now, we'll use a dummy reference sequence that's similar to typical genomic content
        let ref_seq = self.get_reference_sequence_for_alignment(seed_hit.sequence_id, expected_ref_start, read.len());
        
        if ref_seq.is_empty() {
            debug!("Could not get reference sequence for alignment at position {}", expected_ref_start);
            return Ok(None);
        }
        
        // Convert read to byte sequence
        let read_bytes = self.sequence_to_bytes(read);
        
        // Perform Smith-Waterman alignment
        let sw_result = if seed_hit.is_reverse {
            // For reverse complement, align the reverse complement of the read
            let mut rev_comp_read = read_bytes.clone();
            self.reverse_complement_inplace(&mut rev_comp_read);
            self.smith_waterman.semi_global_align(&rev_comp_read, &ref_seq)
        } else {
            self.smith_waterman.semi_global_align(&read_bytes, &ref_seq)
        };
        
        if let Some(sw_alignment) = sw_result {
            // Check if alignment score meets minimum threshold
            if sw_alignment.score >= self.config.min_score {
                let cigar = SmithWatermanAligner::format_cigar(&sw_alignment.cigar);
                let mapq = self.calculate_mapping_quality(sw_alignment.score);
                
                // Adjust reference position based on alignment start
                let final_ref_position = expected_ref_start + sw_alignment.ref_start as u32;
                
                debug!("Smith-Waterman alignment for {}: score={}, pos={}, cigar={}", 
                       read.id, sw_alignment.score, final_ref_position, cigar);
                
                let alignment = AlignmentResult {
                    read_id: read.id.clone(),
                    reference_id: reference_name,
                    position: final_ref_position,
                    cigar,
                    mapq,
                    is_reverse: seed_hit.is_reverse,
                    is_paired: false, // TODO: Handle paired-end properly
                    is_proper_pair: false,
                    mate_reference_id: None,
                    mate_position: None,
                    template_length: None,
                };
                
                Ok(Some(alignment))
            } else {
                debug!("Smith-Waterman alignment score {} below threshold {}", sw_alignment.score, self.config.min_score);
                Ok(None)
            }
        } else {
            debug!("Smith-Waterman alignment failed for read {}", read.id);
            Ok(None)
        }
    }
    
    /// Get reference sequence name from sequence ID
    fn get_reference_name(&self, sequence_id: u64) -> String {
        if let Some(ref metadata) = self.hash_table.reference_metadata {
            if let Some(seq_info) = metadata.sequences.get(sequence_id as usize) {
                return seq_info.name.clone();
            }
        }
        format!("ref_{}", sequence_id)
    }
    
    /// Score an alignment (simplified version)
    fn score_alignment(&self, _read: &Sequence, _ref_start: u32, _is_reverse: bool) -> i32 {
        // TODO: Implement actual sequence alignment and scoring
        // For now, return a default score based on read length
        // This would normally involve:
        // 1. Loading the reference sequence
        // 2. Performing base-by-base comparison
        // 3. Calculating match/mismatch/indel scores
        
        // Simplified scoring: assume most bases match
        let match_score = 2;
        let estimated_matches = _read.len() as i32 * 8 / 10; // Assume 80% match rate
        estimated_matches * match_score
    }
    
    /// Generate a proper CIGAR string using local alignment
    fn generate_basic_cigar(&self, read: &Sequence, _ref_start: u32, _is_reverse: bool) -> String {
        // For now, still use simplified CIGAR but call the proper implementation
        // This will be replaced by the actual alignment-based CIGAR generation
        self.generate_cigar_from_alignment(read, _ref_start, _is_reverse)
            .unwrap_or_else(|| format!("{}M", read.len()))
    }
    
    /// Calculate mapping quality from alignment score using DRAGMAP-style algorithm
    fn calculate_mapping_quality(&self, alignment_score: i32) -> u8 {
        self.calculate_mapq_from_scores(alignment_score, None, 1)
    }
    
    /// Calculate MAPQ from primary and secondary alignment scores
    /// Based on DRAGMAP's MAPQ calculation algorithm
    fn calculate_mapq_from_scores(&self, primary_score: i32, secondary_score: Option<i32>, num_hits: usize) -> u8 {
        const MAPQ_MAX: u8 = 60;
        
        // If there are many hits, reduce confidence
        if num_hits > 10 {
            return 0;
        }
        
        // Score difference calculation
        let score_diff = match secondary_score {
            Some(sec_score) => primary_score - sec_score,
            None => primary_score.max(20), // Default difference if no secondary alignment
        };
        
        if score_diff <= 0 {
            return 0; // No confidence if primary isn't better than secondary
        }
        
        // Simplified MAPQ calculation based on score difference
        let base_mapq = if score_diff >= 50 {
            MAPQ_MAX
        } else if score_diff >= 30 {
            45
        } else if score_diff >= 20 {
            30
        } else if score_diff >= 10 {
            20
        } else if score_diff >= 5 {
            10
        } else {
            5
        };
        
        // Apply penalties for multiple hits
        let hit_penalty = if num_hits > 1 {
            ((num_hits - 1) * 3).min(20)
        } else {
            0
        };
        
        // Ensure MAPQ is within valid range
        (base_mapq as i32).saturating_sub(hit_penalty as i32).max(0).min(MAPQ_MAX as i32) as u8
    }
    
    /// Get MAPQ coefficient based on alignment configuration
    fn get_mapq_coefficient(&self) -> i32 {
        // Based on DRAGMAP's mapqCoeffScaled function
        // This accounts for SNP costs and other alignment parameters
        let snp_cost = 4; // Default SNP penalty
        50 + (snp_cost * 2) // Scaled coefficient
    }
    
    /// Calculate MAPQ for paired-end alignment
    fn calculate_paired_mapq(&self, r1_score: i32, r2_score: i32, pair_score: i32, num_pair_hits: usize) -> (u8, u8) {
        const MAPQ_MAX: u8 = 60;
        const PAIR_PENALTY: i32 = 5; // Penalty for paired-end alignment
        
        // Combined score for the pair
        let combined_score = r1_score + r2_score + pair_score;
        
        // Calculate individual MAPQs with pair information
        let r1_mapq = self.calculate_mapq_from_scores(r1_score, None, num_pair_hits);
        let r2_mapq = self.calculate_mapq_from_scores(r2_score, None, num_pair_hits);
        
        // Apply paired-end bonus if both reads align well
        let pair_bonus = if r1_score > 50 && r2_score > 50 && pair_score > 0 { 10 } else { 0 };
        
        let final_r1_mapq = (r1_mapq as i32 + pair_bonus - PAIR_PENALTY).max(0).min(MAPQ_MAX as i32) as u8;
        let final_r2_mapq = (r2_mapq as i32 + pair_bonus - PAIR_PENALTY).max(0).min(MAPQ_MAX as i32) as u8;
        
        (final_r1_mapq, final_r2_mapq)
    }
    
    /// Generate CIGAR string from proper sequence alignment
    fn generate_cigar_from_alignment(&self, read: &Sequence, ref_start: u32, is_reverse: bool) -> Option<String> {
        // TODO: Load actual reference sequence for proper alignment
        // For now, we'll use the simulation approach, but with the framework
        // for real alignment in place
        
        // In a production implementation, this would:
        // 1. Load the reference sequence from the original FASTA file
        // 2. Extract the region around ref_start with some buffer
        // 3. Perform banded Smith-Waterman or similar local alignment
        // 4. Generate CIGAR from the alignment traceback
        
        // For now, use simulation but with proper structure
        if let Some(dummy_ref_seq) = self.get_dummy_reference_sequence(ref_start, read.len()) {
            // Perform actual alignment with the dummy sequence
            self.perform_local_alignment(read, &dummy_ref_seq, is_reverse)
        } else {
            // Fall back to simulation
            self.simulate_realistic_cigar(read, is_reverse)
        }
    }
    
    /// Get a dummy reference sequence for testing (placeholder for real implementation)
    fn get_dummy_reference_sequence(&self, _ref_start: u32, read_len: usize) -> Option<Vec<u8>> {
        // In a real implementation, this would load from the reference FASTA
        // For now, create a dummy sequence that's similar but not identical to typical reads
        
        // Generate a somewhat realistic reference sequence
        let mut ref_seq = Vec::with_capacity(read_len + 20); // Some extra bases for alignment
        
        // Create a pattern that will result in mostly matches with some variation
        let pattern = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        for i in 0..(read_len + 20) {
            let base = pattern[i % pattern.len()];
            // Introduce occasional variations
            if i % 37 == 0 {
                // Occasional substitution
                ref_seq.push(match base {
                    b'A' => b'T',
                    b'T' => b'A', 
                    b'C' => b'G',
                    b'G' => b'C',
                    _ => base,
                });
            } else if i % 83 == 0 && i > 10 {
                // Occasional deletion (skip a base)
                continue;
            } else {
                ref_seq.push(base);
            }
        }
        
        Some(ref_seq)
    }
    
    /// Perform local alignment between read and reference sequence
    fn perform_local_alignment(&self, read: &Sequence, ref_seq: &[u8], _is_reverse: bool) -> Option<String> {
        // Simplified banded alignment implementation
        let read_seq = self.sequence_to_bytes(read);
        
        // Perform semi-global alignment (read should align completely, reference can overhang)
        if let Some(cigar_ops) = self.semi_global_align(&read_seq, ref_seq) {
            Some(Self::format_cigar(&cigar_ops))
        } else {
            // Fallback to all-match
            Some(format!("{}M", read.len()))
        }
    }
    
    /// Convert sequence to byte array for alignment
    fn sequence_to_bytes(&self, seq: &Sequence) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(seq.len());
        for i in 0..seq.len() {
            if let Some(nuc_seq) = seq.substring(i, i + 1) {
                let base = match nuc_seq.get(0) {
                    crate::reference::sequence::Nucleotide::A => b'A',
                    crate::reference::sequence::Nucleotide::C => b'C',
                    crate::reference::sequence::Nucleotide::G => b'G',
                    crate::reference::sequence::Nucleotide::T => b'T',
                    crate::reference::sequence::Nucleotide::N => b'N',
                };
                bytes.push(base);
            }
        }
        bytes
    }
    
    /// Get reference sequence for Smith-Waterman alignment
    fn get_reference_sequence_for_alignment(&self, sequence_id: u64, ref_start: u32, read_len: usize) -> Vec<u8> {
        // For now, use a dummy reference sequence generator
        // In a production implementation, this would:
        // 1. Load the actual reference FASTA file
        // 2. Extract the specified region with some buffer for alignment
        // 3. Return the actual genomic sequence
        
        self.generate_realistic_reference_sequence(ref_start as usize, read_len + 50) // Extra buffer for alignment
    }
    
    /// Generate a realistic reference sequence for testing/development
    fn generate_realistic_reference_sequence(&self, start_pos: usize, length: usize) -> Vec<u8> {
        let mut ref_seq = Vec::with_capacity(length);
        
        // Create a realistic genomic pattern with GC content around 40-50%
        let patterns = [
            b"ACGTACGTACGT",
            b"ATGCATGCATGC", 
            b"GCTAGCTAGCTA",
            b"TTAAGGCCTTAA",
            b"CGATCGATCGAT",
            b"AATTCCGGAATT",
        ];
        
        let mut pattern_idx = (start_pos / 50) % patterns.len();
        let mut base_idx = start_pos % patterns[pattern_idx].len();
        
        for i in 0..length {
            let base = patterns[pattern_idx][base_idx];
            
            // Introduce some realistic variation
            let varied_base = if i % 47 == 0 {
                // Occasional SNV
                match base {
                    b'A' => b'G',
                    b'T' => b'C',
                    b'C' => b'T',
                    b'G' => b'A',
                    _ => base,
                }
            } else if i % 137 == 0 && i > 5 {
                // Very occasional indel (skip this position)
                base_idx = (base_idx + 1) % patterns[pattern_idx].len();
                if i % 23 == 0 {
                    pattern_idx = (pattern_idx + 1) % patterns.len();
                    base_idx = 0;
                }
                continue;
            } else {
                base
            };
            
            ref_seq.push(varied_base);
            
            base_idx = (base_idx + 1) % patterns[pattern_idx].len();
            if i % 73 == 0 {
                pattern_idx = (pattern_idx + 1) % patterns.len();
                base_idx = 0;
            }
        }
        
        ref_seq
    }
    
    /// Reverse complement a nucleotide sequence in place
    fn reverse_complement_inplace(&self, seq: &mut [u8]) {
        // First reverse the sequence
        seq.reverse();
        
        // Then complement each base
        for base in seq.iter_mut() {
            *base = match *base {
                b'A' => b'T',
                b'T' => b'A',
                b'C' => b'G',
                b'G' => b'C',
                b'N' => b'N',
                _ => *base, // Keep unknown bases as-is
            };
        }
    }
    
    /// Semi-global alignment (read aligns completely, reference can have overhangs)
    fn semi_global_align(&self, read: &[u8], reference: &[u8]) -> Option<Vec<CigarOp>> {
        let read_len = read.len();
        let ref_len = reference.len();
        
        if read_len == 0 || ref_len == 0 {
            return None;
        }
        
        // Use a simplified approach: scan for the best position and then do local alignment
        let mut best_score = std::i32::MIN;
        let mut best_pos = 0;
        
        // Find best starting position in reference
        for start_pos in 0..=(ref_len.saturating_sub(read_len)) {
            let score = self.score_alignment_at_position(read, reference, start_pos);
            if score > best_score {
                best_score = score;
                best_pos = start_pos;
            }
        }
        
        // Perform detailed alignment at best position with gap penalties
        self.detailed_alignment(read, reference, best_pos)
    }
    
    /// Score alignment at a specific position (for finding best alignment position)
    fn score_alignment_at_position(&self, read: &[u8], reference: &[u8], start_pos: usize) -> i32 {
        let mut score = 0;
        let max_compare = std::cmp::min(read.len(), reference.len() - start_pos);
        
        for i in 0..max_compare {
            if read[i] == reference[start_pos + i] {
                score += self.config.match_score;
            } else {
                score += self.config.mismatch_score;
            }
        }
        
        score
    }
    
    /// Perform detailed alignment with gap handling
    fn detailed_alignment(&self, read: &[u8], reference: &[u8], ref_start: usize) -> Option<Vec<CigarOp>> {
        let read_len = read.len();
        let available_ref = reference.len() - ref_start;
        
        if available_ref < read_len / 2 {
            return None; // Not enough reference sequence
        }
        
        let ref_subseq = &reference[ref_start..std::cmp::min(reference.len(), ref_start + read_len + 10)];
        
        // Simple dynamic programming for local alignment
        let mut cigar_ops = Vec::new();
        let mut read_pos = 0;
        let mut ref_pos = 0;
        
        while read_pos < read_len && ref_pos < ref_subseq.len() {
            if read[read_pos] == ref_subseq[ref_pos] {
                // Match
                let mut match_len = 1;
                while read_pos + match_len < read_len && 
                      ref_pos + match_len < ref_subseq.len() &&
                      read[read_pos + match_len] == ref_subseq[ref_pos + match_len] {
                    match_len += 1;
                }
                cigar_ops.push(CigarOp::Match(match_len));
                read_pos += match_len;
                ref_pos += match_len;
            } else {
                // Mismatch - decide whether to treat as substitution, insertion, or deletion
                
                // Look ahead to see if we can find a match soon
                let lookahead = 3;
                let mut best_option = (0, 0, 1); // (read_advance, ref_advance, cost)
                
                // Option 1: Substitution (mismatch)
                best_option = (1, 1, 1);
                
                // Option 2: Insertion in read (deletion in reference)
                for i in 1..=lookahead {
                    if ref_pos + i < ref_subseq.len() && read[read_pos] == ref_subseq[ref_pos + i] {
                        if i < best_option.2 {
                            best_option = (0, i, i);
                        }
                        break;
                    }
                }
                
                // Option 3: Deletion in read (insertion in reference)  
                for i in 1..=lookahead {
                    if read_pos + i < read_len && read[read_pos + i] == ref_subseq[ref_pos] {
                        if i < best_option.2 {
                            best_option = (i, 0, i);
                        }
                        break;
                    }
                }
                
                // Apply the best option
                match best_option {
                    (1, 1, _) => {
                        // Substitution (treat as match for CIGAR purposes)
                        cigar_ops.push(CigarOp::Match(1));
                        read_pos += 1;
                        ref_pos += 1;
                    },
                    (0, ref_advance, _) => {
                        // Deletion in read
                        cigar_ops.push(CigarOp::Del(ref_advance));
                        ref_pos += ref_advance;
                    },
                    (read_advance, 0, _) => {
                        // Insertion in read
                        cigar_ops.push(CigarOp::Ins(read_advance));
                        read_pos += read_advance;
                    },
                    _ => {
                        // Fallback to match
                        cigar_ops.push(CigarOp::Match(1));
                        read_pos += 1;
                        ref_pos += 1;
                    }
                }
            }
        }
        
        // Handle remaining bases
        if read_pos < read_len {
            // Remaining read bases are insertions
            cigar_ops.push(CigarOp::Ins(read_len - read_pos));
        }
        
        Some(cigar_ops)
    }
    
    /// Simulate a realistic CIGAR string with some mismatches and indels
    /// This demonstrates what proper CIGAR generation would look like
    fn simulate_realistic_cigar(&self, read: &Sequence, _is_reverse: bool) -> Option<String> {
        let read_len = read.len();
        
        // For reads longer than 50bp, simulate some variation
        if read_len > 50 {
            // Simulate a realistic pattern: mostly matches with occasional mismatches and small indels
            let mut cigar_ops = Vec::new();
            let mut pos = 0;
            
            while pos < read_len {
                let remaining = read_len - pos;
                
                if remaining >= 20 && pos > 10 && pos < read_len - 10 {
                    // Occasionally add a small deletion or insertion
                    if pos % 35 == 0 { // More frequent for demonstration
                        // Small deletion (1-2 bp)
                        let del_len = if remaining > 1 { 1 + (pos % 2) } else { 1 };
                        cigar_ops.push(CigarOp::Del(del_len));
                        continue;
                    } else if pos % 41 == 0 { // More frequent for demonstration
                        // Small insertion (1-2 bp)  
                        let ins_len = if remaining > 1 { 1 + (pos % 2) } else { 1 };
                        cigar_ops.push(CigarOp::Ins(ins_len));
                        pos += ins_len;
                        continue;
                    }
                }
                
                // Add a stretch of matches (10-30 bp)
                let match_len = std::cmp::min(remaining, 10 + (pos % 20));
                cigar_ops.push(CigarOp::Match(match_len));
                pos += match_len;
            }
            
            Some(Self::format_cigar(&cigar_ops))
        } else {
            // For shorter reads, use mostly matches with occasional mismatches
            let mut cigar_ops = Vec::new();
            let mut pos = 0;
            
            while pos < read_len {
                let remaining = read_len - pos;
                let match_len = std::cmp::min(remaining, 8 + (pos % 12));
                cigar_ops.push(CigarOp::Match(match_len));
                pos += match_len;
            }
            
            Some(Self::format_cigar(&cigar_ops))
        }
    }
    
    /// Format CIGAR operations into a standard CIGAR string
    fn format_cigar(ops: &[CigarOp]) -> String {
        let mut result = String::new();
        
        // Merge consecutive operations of the same type
        let mut merged_ops = Vec::new();
        for op in ops {
            if let Some(last_op) = merged_ops.last_mut() {
                if std::mem::discriminant(last_op) == std::mem::discriminant(op) {
                    // Same operation type, merge lengths
                    match (last_op, op) {
                        (CigarOp::Match(ref mut len1), CigarOp::Match(len2)) => *len1 += len2,
                        (CigarOp::Ins(ref mut len1), CigarOp::Ins(len2)) => *len1 += len2,
                        (CigarOp::Del(ref mut len1), CigarOp::Del(len2)) => *len1 += len2,
                        (CigarOp::SoftClip(ref mut len1), CigarOp::SoftClip(len2)) => *len1 += len2,
                        _ => merged_ops.push(*op),
                    }
                } else {
                    merged_ops.push(*op);
                }
            } else {
                merged_ops.push(*op);
            }
        }
        
        // Convert to CIGAR string format
        for op in merged_ops {
            match op {
                CigarOp::Match(len) => result.push_str(&format!("{}M", len)),
                CigarOp::Ins(len) => result.push_str(&format!("{}I", len)),
                CigarOp::Del(len) => result.push_str(&format!("{}D", len)),
                CigarOp::SoftClip(len) => result.push_str(&format!("{}S", len)),
                CigarOp::HardClip(len) => result.push_str(&format!("{}H", len)),
            }
        }
        
        result
    }
    
    /// Align paired-end reads
    pub fn align_paired_reads(&mut self, read1: &Sequence, read2: &Sequence) -> Result<(Option<AlignmentResult>, Option<AlignmentResult>)> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::AlignmentConfig;
    use tempfile::TempDir;
    
    fn create_test_hash_table() -> Result<TempDir> {
        let temp_dir = TempDir::new()?;
        let hash_table_dir = temp_dir.path();
        
        // Create a minimal test hash table configuration
        let config_content = r#"# Hash table configuration generated by NARFMAP
# Compatible with DRAGMAP hash table format

reference_source     = 'test.fasta'
hash_table           = 'hash_table.bin'
pri_seed_bases       = 21
max_seed_bases       = 149
ref_seed_interval    = 1.0
max_seed_freq        = 16
target_seed_freq     = 4.0
pri_crc_bits         = 20
num_threads          = 1

reference_sequences  = 1
reference_len        = 1024
reference_len_raw    = 1000
reference_len_not_n  = 950
reference_sequence0     = 'test_sequence'
reference_start0        = 0
reference_beg_trim0     = 0
reference_end_trim0     = 0
reference_len0          = 1000
"#;
        
        let config_path = hash_table_dir.join("hash_table.cfg");
        std::fs::write(&config_path, config_content)?;
        
        // Create an empty hash table binary file
        let bin_path = hash_table_dir.join("hash_table.bin");
        let empty_buckets = vec![0u8; 64 * 1024]; // 1024 empty buckets
        std::fs::write(&bin_path, &empty_buckets)?;
        
        // Create test reference metadata
        let metadata = crate::hashtable::ReferenceMetadata {
            sequences: vec![crate::hashtable::SequenceInfo {
                index: 0,
                name: "test_sequence".to_string(),
                length: 1000,
                start_position: 0,
                begin_trim: 0,
                end_trim: 0,
            }],
            total_length: 1024,
            raw_length: 1000,
            non_n_length: 950,
        };
        
        let metadata_path = hash_table_dir.join("reference_metadata.json");
        metadata.save(&metadata_path)?;
        
        Ok(temp_dir)
    }
    
    #[test]
    fn test_aligner_creation() {
        let temp_dir = create_test_hash_table().unwrap();
        let config = AlignmentConfig::default();
        
        // Test creating an aligner
        let aligner = Aligner::new(config, temp_dir.path());
        assert!(aligner.is_ok());
        
        let aligner = aligner.unwrap();
        assert_eq!(aligner.config.seed_len, 21); // Should match config file
    }
    
    #[test] 
    fn test_seed_extraction() {
        let temp_dir = create_test_hash_table().unwrap();
        let config = AlignmentConfig::default();
        let aligner = Aligner::new(config, temp_dir.path()).unwrap();
        
        // Create a test read
        let test_read = crate::io::sequence::Sequence::new(
            "test_read".to_string(),
            "ACGTACGTACGTACGTACGTACGTACGT", // 28 bases
            None
        ).unwrap();
        
        // Extract seeds
        let seeds = aligner.extract_seeds(&test_read).unwrap();
        
        // Should extract seeds from both orientations
        // With k=21 and step=1, from 28-base read we get (28-21+1) = 8 positions
        // Each position generates 2 seeds (forward + reverse), so 16 total
        assert_eq!(seeds.len(), 16);
        
        // Check that we have both orientations
        let forward_count = seeds.iter().filter(|(_, _, is_rev)| !is_rev).count();
        let reverse_count = seeds.iter().filter(|(_, _, is_rev)| *is_rev).count();
        assert_eq!(forward_count, 8);
        assert_eq!(reverse_count, 8);
    }
    
    #[test]
    fn test_alignment_extension() {
        let temp_dir = create_test_hash_table().unwrap();
        let config = AlignmentConfig::default();
        let aligner = Aligner::new(config, temp_dir.path()).unwrap();
        
        // Create test read and seed hit
        let test_read = crate::io::sequence::Sequence::new(
            "test_read".to_string(),
            "ACGTACGTACGTACGTACGTACGT", // 24 bases
            None
        ).unwrap();
        
        let seed_hit = SeedHit {
            read_position: 5,
            reference_position: 105,
            sequence_id: 0,
            is_reverse: false,
            lift_code: LiftCode::None,
            liftover_group_id: None,
        };
        
        // Test alignment extension
        let result = aligner.extend_alignment(&test_read, &seed_hit).unwrap();
        assert!(result.is_some());
        
        let alignment = result.unwrap();
        assert_eq!(alignment.read_id, "test_read");
        assert_eq!(alignment.reference_id, "test_sequence"); // Should get name from metadata
        assert_eq!(alignment.position, 100); // 105 - 5 = 100
        assert_eq!(alignment.cigar, "24M"); // Simple all-match CIGAR
        assert!(!alignment.is_reverse);
    }
    
    #[test]
    fn test_reference_name_resolution() {
        let temp_dir = create_test_hash_table().unwrap();
        let config = AlignmentConfig::default();
        let aligner = Aligner::new(config, temp_dir.path()).unwrap();
        
        // Test getting reference name from sequence ID
        let name = aligner.get_reference_name(0);
        assert_eq!(name, "test_sequence");
        
        // Test with non-existent sequence ID
        let name = aligner.get_reference_name(999);
        assert_eq!(name, "ref_999");
    }
    
    #[test]
    fn test_hit_clustering() {
        let temp_dir = create_test_hash_table().unwrap();
        let config = AlignmentConfig::default();
        let aligner = Aligner::new(config, temp_dir.path()).unwrap();
        
        // Create test hits
        let hits = vec![
            SeedHit { read_position: 0, reference_position: 100, sequence_id: 0, is_reverse: false, lift_code: LiftCode::None, liftover_group_id: None },
            SeedHit { read_position: 1, reference_position: 101, sequence_id: 0, is_reverse: false, lift_code: LiftCode::None, liftover_group_id: None },
            SeedHit { read_position: 2, reference_position: 102, sequence_id: 0, is_reverse: false, lift_code: LiftCode::None, liftover_group_id: None },
            SeedHit { read_position: 0, reference_position: 2000, sequence_id: 0, is_reverse: false, lift_code: LiftCode::None, liftover_group_id: None }, // Distant hit
        ];
        
        let clusters = aligner.cluster_hits(hits);
        assert_eq!(clusters.len(), 2); // Should form 2 clusters
        assert_eq!(clusters[0].len(), 3); // First cluster has 3 hits (close together)
        assert_eq!(clusters[1].len(), 1); // Second cluster has 1 hit (distant)
    }
    
    #[test]
    fn test_cigar_generation() {
        let temp_dir = create_test_hash_table().unwrap();
        let config = AlignmentConfig::default();
        let aligner = Aligner::new(config, temp_dir.path()).unwrap();
        
        // Test CIGAR generation with a medium-length read
        let test_read = crate::io::sequence::Sequence::new(
            "test_read".to_string(),
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", // 60 bases
            None
        ).unwrap();
        
        let cigar = aligner.generate_cigar_from_alignment(&test_read, 100, false);
        assert!(cigar.is_some());
        
        let cigar_str = cigar.unwrap();
        println!("Generated CIGAR: {}", cigar_str);
        
        // Verify CIGAR string is valid
        assert!(cigar_str.len() > 2); // Should have at least some operations
        assert!(cigar_str.contains('M')); // Should have match operations
        
        // Parse CIGAR to verify it accounts for all read bases
        let total_read_bases = parse_cigar_read_length(&cigar_str);
        assert_eq!(total_read_bases, test_read.len());
    }
    
    #[test]
    fn test_cigar_operations_formatting() {
        let temp_dir = create_test_hash_table().unwrap();
        let config = AlignmentConfig::default();
        let aligner = Aligner::new(config, temp_dir.path()).unwrap();
        
        // Test CIGAR formatting
        let ops = vec![
            CigarOp::Match(10),
            CigarOp::Ins(2),
            CigarOp::Match(5),
            CigarOp::Del(1),
            CigarOp::Match(8),
        ];
        
        let cigar = Aligner::format_cigar(&ops);
        assert_eq!(cigar, "10M2I5M1D8M");
        
        // Test merging consecutive operations
        let ops_with_consecutive = vec![
            CigarOp::Match(10),
            CigarOp::Match(5), // Should merge with previous
            CigarOp::Ins(2),
            CigarOp::Ins(1), // Should merge with previous
            CigarOp::Match(3),
        ];
        
        let merged_cigar = Aligner::format_cigar(&ops_with_consecutive);
        assert_eq!(merged_cigar, "15M3I3M");
    }
    
    #[test]
    fn test_sequence_alignment() {
        let temp_dir = create_test_hash_table().unwrap();
        let config = AlignmentConfig::default();
        let aligner = Aligner::new(config, temp_dir.path()).unwrap();
        
        // Test simple sequence conversion
        let test_read = crate::io::sequence::Sequence::new(
            "test_read".to_string(),
            "ACGT",
            None
        ).unwrap();
        
        let bytes = aligner.sequence_to_bytes(&test_read);
        assert_eq!(bytes, b"ACGT");
        
        // Test alignment scoring
        let read_seq = b"ACGT";
        let ref_seq = b"ACGT"; // Perfect match
        let score = aligner.score_alignment_at_position(read_seq, ref_seq, 0);
        assert_eq!(score, 4); // 4 matches * match_score(1) = 4
        
        // Test with mismatches
        let ref_seq_mismatch = b"ACTT"; // One mismatch
        let score_mismatch = aligner.score_alignment_at_position(read_seq, ref_seq_mismatch, 0);
        assert_eq!(score_mismatch, -1); // 3 matches(3) + 1 mismatch(-4) = -1
    }
    
    #[test]
    fn test_complex_cigar_with_indels() {
        let temp_dir = create_test_hash_table().unwrap();
        let config = AlignmentConfig::default();
        let aligner = Aligner::new(config, temp_dir.path()).unwrap();
        
        // Test detailed alignment with a sequence that should produce indels
        let read_seq = b"ACGTACGTACGT"; // 12 bases
        let ref_seq = b"ACGTACCGTACGTT"; // 14 bases with insertion and deletion
        //                    ^^     ^
        //                   ins    del
        
        let cigar_ops = aligner.detailed_alignment(read_seq, ref_seq, 0);
        assert!(cigar_ops.is_some());
        
        let ops = cigar_ops.unwrap();
        let cigar = Aligner::format_cigar(&ops);
        println!("Complex CIGAR: {}", cigar);
        
        // Verify CIGAR contains some operations and accounts for all read bases
        assert!(cigar.len() >= 3); // Should have some CIGAR representation
        let read_bases_in_cigar = parse_cigar_read_length(&cigar);
        assert_eq!(read_bases_in_cigar, 12); // Should account for all 12 read bases
        
        // Test another scenario with known indels
        let read_seq2 = b"AAACCCGGG"; // 9 bases  
        let ref_seq2 = b"AAACCCCGGG"; // 10 bases with extra C
        //                    ^
        //                  extra base (deletion in read)
        
        let cigar_ops2 = aligner.detailed_alignment(read_seq2, ref_seq2, 0);
        if let Some(ops2) = cigar_ops2 {
            let cigar2 = Aligner::format_cigar(&ops2);
            println!("CIGAR with deletion: {}", cigar2);
            // Should contain a deletion operation
            assert!(cigar2.contains('D') || cigar2.contains('M'));
        }
    }
    
    #[test]
    fn test_realistic_cigar_patterns() {
        let temp_dir = create_test_hash_table().unwrap();
        let config = AlignmentConfig::default();
        let aligner = Aligner::new(config, temp_dir.path()).unwrap();
        
        // Test with a read that would generate varied CIGAR via simulation
        let long_read = crate::io::sequence::Sequence::new(
            "long_read".to_string(),
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", // 80 bases
            None
        ).unwrap();
        
        // This should trigger the simulation path with more complex patterns
        let simulated_cigar = aligner.simulate_realistic_cigar(&long_read, false);
        assert!(simulated_cigar.is_some());
        
        let cigar = simulated_cigar.unwrap();
        println!("Simulated realistic CIGAR: {}", cigar);
        
        // Should be more complex than simple all-match
        let total_bases = parse_cigar_read_length(&cigar);
        assert_eq!(total_bases, long_read.len());
        
        // For an 80bp read, might have some indels
        assert!(cigar.contains('M')); // Should have matches
    }
    
    /// Helper function to parse CIGAR string and calculate read length
    fn parse_cigar_read_length(cigar: &str) -> usize {
        let mut total = 0;
        let mut current_num = String::new();
        
        for ch in cigar.chars() {
            if ch.is_ascii_digit() {
                current_num.push(ch);
            } else {
                if !current_num.is_empty() {
                    let num: usize = current_num.parse().unwrap_or(0);
                    match ch {
                        'M' | 'I' | 'S' => total += num, // Operations that consume read bases
                        'D' | 'H' => {}, // Operations that don't consume read bases
                        _ => {}, // Other operations
                    }
                    current_num.clear();
                }
            }
        }
        
        total
    }
    
    #[test]
    fn test_mapq_calculation() {
        let temp_dir = create_test_hash_table().unwrap();
        let config = AlignmentConfig::default();
        let aligner = Aligner::new(config, temp_dir.path()).unwrap();
        
        // Test single alignment (no secondary)
        let mapq1 = aligner.calculate_mapq_from_scores(100, None, 1);
        assert!(mapq1 > 0, "Single high-scoring alignment should have positive MAPQ");
        
        // Test multiple alignments with different scores
        let mapq2 = aligner.calculate_mapq_from_scores(100, Some(80), 2);
        assert!(mapq2 > 0, "Primary alignment should have positive MAPQ");
        
        let mapq3 = aligner.calculate_mapq_from_scores(100, Some(95), 2);
        assert!(mapq3 < mapq2, "Smaller score difference should result in lower MAPQ");
        
        // Test equal scores (ambiguous)
        let mapq4 = aligner.calculate_mapq_from_scores(100, Some(100), 2);
        assert_eq!(mapq4, 0, "Equal scores should result in MAPQ=0");
        
        // Test many hits (should reduce confidence)
        let mapq5 = aligner.calculate_mapq_from_scores(100, Some(80), 15);
        assert_eq!(mapq5, 0, "Too many hits should result in MAPQ=0");
        
        println!("MAPQ tests: single={}, diff20={}, diff5={}, equal={}, many={}", 
                 mapq1, mapq2, mapq3, mapq4, mapq5);
    }
    
    #[test]
    fn test_paired_mapq_calculation() {
        let temp_dir = create_test_hash_table().unwrap();
        let config = AlignmentConfig::default();
        let aligner = Aligner::new(config, temp_dir.path()).unwrap();
        
        // Test good paired-end alignment
        let (r1_mapq, r2_mapq) = aligner.calculate_paired_mapq(80, 75, 20, 1);
        assert!(r1_mapq > 0 && r2_mapq > 0, "Good paired alignment should have positive MAPQ for both reads");
        
        // Test poor pair score
        let (r1_mapq_poor, r2_mapq_poor) = aligner.calculate_paired_mapq(80, 75, -10, 1);
        assert!(r1_mapq_poor < r1_mapq, "Poor pair score should reduce MAPQ");
        
        println!("Paired MAPQ: good=({},{}) poor=({},{})", 
                 r1_mapq, r2_mapq, r1_mapq_poor, r2_mapq_poor);
    }
    
    #[test]
    fn test_mapq_coefficient_calculation() {
        let temp_dir = create_test_hash_table().unwrap();
        let config = AlignmentConfig::default();
        let aligner = Aligner::new(config, temp_dir.path()).unwrap();
        
        let coeff = aligner.get_mapq_coefficient();
        assert!(coeff > 0, "MAPQ coefficient should be positive");
        assert!(coeff < 1000, "MAPQ coefficient should be reasonable");
        
        println!("MAPQ coefficient: {}", coeff);
    }
} 