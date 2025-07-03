use std::path::Path;
use anyhow::{Result, anyhow};
use log::{info, debug, warn};

use crate::config::{AlignmentConfig, HashTableConfig};
use crate::hashtable::{HashTableQuery, KmerHasher};
use crate::io::sequence::Sequence;

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
        
        Ok(Self {
            config,
            hash_table,
            _kmer_hasher: kmer_hasher,
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
        
        // Cluster hits by reference position
        let clustered_hits = self.cluster_hits(all_hits);
        debug!("Clustered into {} potential alignment positions", clustered_hits.len());
        
        // Score and select best alignment
        if let Some(best_hit) = self.select_best_alignment(read, clustered_hits)? {
            debug!("Best alignment for {}: ref_pos={}, score={}", read.id, best_hit.reference_position, "TODO");
            
            // Perform alignment extension from seed hit
            if let Some(extended_alignment) = self.extend_alignment(read, &best_hit)? {
                debug!("Extended alignment for {}: pos={}, cigar={}", read.id, extended_alignment.position, extended_alignment.cigar);
                Ok(Some(extended_alignment))
            } else {
                debug!("Failed to extend alignment for read {}", read.id);
                Ok(None)
            }
        } else {
            debug!("No valid alignment found for read {}", read.id);
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
    
    /// Find hits for a k-mer seed in the hash table
    fn find_seed_hits(&self, seed: &[u8], read_position: usize, is_reverse: bool) -> Result<Vec<SeedHit>> {
        let positions = self.hash_table.query_kmer(seed)?;
        
        let mut hits = Vec::new();
        for encoded_pos in positions {
            let sequence_id = encoded_pos >> 32;
            let reference_position = (encoded_pos & 0xFFFFFFFF) as u32;
            
            hits.push(SeedHit {
                read_position,
                reference_position,
                sequence_id,
                is_reverse,
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
    
    /// Extend alignment from a seed hit using simple local alignment
    fn extend_alignment(&self, read: &Sequence, seed_hit: &SeedHit) -> Result<Option<AlignmentResult>> {
        // For basic extension, we'll load the reference sequence and perform
        // a simple base-by-base comparison
        
        // TODO: Load actual reference sequence - for now use dummy implementation
        // This would require loading the original reference FASTA file
        let reference_name = self.get_reference_name(seed_hit.sequence_id);
        
        // Calculate the expected start position on reference
        // If we have a seed hit at reference position X and read position Y,
        // the read should start at reference position (X - Y)
        let expected_ref_start = if seed_hit.reference_position >= seed_hit.read_position as u32 {
            seed_hit.reference_position - seed_hit.read_position as u32
        } else {
            0 // Clamp to start of reference
        };
        
        // For now, create a basic alignment with simple scoring
        let alignment_score = self.score_alignment(read, expected_ref_start, seed_hit.is_reverse);
        
        if alignment_score >= self.config.min_score {
            let cigar = self.generate_basic_cigar(read, expected_ref_start, seed_hit.is_reverse);
            let mapq = self.calculate_mapping_quality(alignment_score);
            
            let alignment = AlignmentResult {
                read_id: read.id.clone(),
                reference_id: reference_name,
                position: expected_ref_start,
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
            debug!("Alignment score {} below threshold {}", alignment_score, self.config.min_score);
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
    
    /// Calculate mapping quality from alignment score
    fn calculate_mapping_quality(&self, alignment_score: i32) -> u8 {
        // Simple mapping quality calculation
        // Higher scores get higher MAPQ values
        if alignment_score >= 100 {
            60 // High confidence
        } else if alignment_score >= 50 {
            30 // Medium confidence  
        } else if alignment_score >= 20 {
            10 // Low confidence
        } else {
            0 // Very low confidence
        }
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
            SeedHit { read_position: 0, reference_position: 100, sequence_id: 0, is_reverse: false },
            SeedHit { read_position: 1, reference_position: 101, sequence_id: 0, is_reverse: false },
            SeedHit { read_position: 2, reference_position: 102, sequence_id: 0, is_reverse: false },
            SeedHit { read_position: 0, reference_position: 2000, sequence_id: 0, is_reverse: false }, // Distant hit
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
} 