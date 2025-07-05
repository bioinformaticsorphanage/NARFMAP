/// Paired-end alignment support with mate rescue scanning
/// 
/// This module implements DRAGMAP-compatible paired-end read alignment
/// including mate rescue scanning within insert size intervals.

use std::collections::HashMap;
use anyhow::Result;
use log::{debug, warn};

use crate::io::sequence::Sequence;
use super::{AlignmentResult, SeedHit, Aligner};

/// Paired-end alignment configuration
#[derive(Debug, Clone)]
pub struct PairedEndConfig {
    /// Minimum insert size for valid pairs
    pub min_insert_size: u32,
    /// Maximum insert size for valid pairs
    pub max_insert_size: u32,
    /// Standard deviation for insert size distribution
    pub insert_size_std: f64,
    /// Mean insert size (updated empirically)
    pub mean_insert_size: f64,
    /// Maximum rescue search distance from mate position
    pub max_rescue_distance: u32,
    /// Minimum MAPQ for anchor read in rescue scanning
    pub min_rescue_anchor_mapq: u8,
    /// Enable mate rescue scanning
    pub enable_mate_rescue: bool,
}

impl Default for PairedEndConfig {
    fn default() -> Self {
        Self {
            min_insert_size: 50,
            max_insert_size: 1000,
            insert_size_std: 100.0,
            mean_insert_size: 500.0,
            max_rescue_distance: 1500,
            min_rescue_anchor_mapq: 20,
            enable_mate_rescue: true,
        }
    }
}

/// Paired-end read alignment result
#[derive(Debug, Clone)]
pub struct PairedAlignment {
    /// First read alignment (R1)
    pub read1: Option<AlignmentResult>,
    /// Second read alignment (R2)
    pub read2: Option<AlignmentResult>,
    /// Insert size (if both reads aligned)
    pub insert_size: Option<i32>,
    /// Whether this is a proper pair
    pub is_proper_pair: bool,
    /// Whether mate rescue was used
    pub rescued: bool,
}

/// Insert size statistics for empirical modeling
#[derive(Debug, Clone)]
pub struct InsertSizeStats {
    /// Observed insert sizes
    sizes: Vec<i32>,
    /// Current mean
    pub mean: f64,
    /// Current standard deviation
    pub std_dev: f64,
    /// Number of observations
    pub count: usize,
    /// Update frequency (recalculate stats every N observations)
    update_frequency: usize,
}

impl InsertSizeStats {
    pub fn new() -> Self {
        Self {
            sizes: Vec::new(),
            mean: 500.0,
            std_dev: 100.0,
            count: 0,
            update_frequency: 1000,
        }
    }
    
    /// Add a new insert size observation
    pub fn observe(&mut self, insert_size: i32) {
        if insert_size > 0 && insert_size < 10000 { // Reasonable bounds
            self.sizes.push(insert_size);
            self.count += 1;
            
            // Update statistics periodically
            if self.count % self.update_frequency == 0 {
                self.update_statistics();
            }
        }
    }
    
    /// Update mean and standard deviation from observations
    fn update_statistics(&mut self) {
        if self.sizes.is_empty() {
            return;
        }
        
        // Calculate mean
        let sum: i32 = self.sizes.iter().sum();
        self.mean = sum as f64 / self.sizes.len() as f64;
        
        // Calculate standard deviation
        let variance: f64 = self.sizes
            .iter()
            .map(|&x| {
                let diff = x as f64 - self.mean;
                diff * diff
            })
            .sum::<f64>() / self.sizes.len() as f64;
        
        self.std_dev = variance.sqrt();
        
        debug!("Updated insert size stats: mean={:.1}, std_dev={:.1}, n={}", 
               self.mean, self.std_dev, self.sizes.len());
        
        // Keep only recent observations to adapt to changing conditions
        if self.sizes.len() > 10000 {
            self.sizes.drain(0..5000); // Keep most recent 5000 observations
        }
    }
    
    /// Check if an insert size is within expected range
    pub fn is_valid_insert_size(&self, insert_size: i32) -> bool {
        let insert_size_f = insert_size as f64;
        let lower_bound = self.mean - 4.0 * self.std_dev;
        let upper_bound = self.mean + 4.0 * self.std_dev;
        insert_size_f >= lower_bound && insert_size_f <= upper_bound
    }
}

/// Paired-end aligner with mate rescue scanning
pub struct PairedEndAligner {
    config: PairedEndConfig,
    insert_stats: InsertSizeStats,
}

impl PairedEndAligner {
    pub fn new(config: PairedEndConfig) -> Self {
        Self {
            config,
            insert_stats: InsertSizeStats::new(),
        }
    }
    
    /// Align a pair of reads with mate rescue scanning
    pub fn align_pair(
        &mut self,
        aligner: &mut Aligner,
        read1: &Sequence,
        read2: &Sequence,
    ) -> Result<PairedAlignment> {
        debug!("Aligning paired reads: {} and {}", read1.id, read2.id);
        
        // Align both reads independently first
        let mut alignment1 = aligner.align_read(read1)?;
        let mut alignment2 = aligner.align_read(read2)?;
        
        // Try to determine proper pairing
        let mut is_proper_pair = false;
        let mut insert_size = None;
        let mut rescued = false;
        
        if let (Some(ref a1), Some(ref a2)) = (&alignment1, &alignment2) {
            // Both reads aligned - check if they form a proper pair
            insert_size = self.calculate_insert_size(a1, a2);
            is_proper_pair = self.is_proper_pair(a1, a2, insert_size);
            
            if let Some(isize) = insert_size {
                self.insert_stats.observe(isize.abs());
            }
        } else if self.config.enable_mate_rescue {
            // One or both reads failed to align - try mate rescue
            let rescue_result = self.attempt_mate_rescue(aligner, read1, read2, &alignment1, &alignment2)?;
            
            if let Some((rescued_a1, rescued_a2)) = rescue_result {
                alignment1 = rescued_a1;
                alignment2 = rescued_a2;
                rescued = true;
                
                if let (Some(ref a1), Some(ref a2)) = (&alignment1, &alignment2) {
                    insert_size = self.calculate_insert_size(a1, a2);
                    is_proper_pair = self.is_proper_pair(a1, a2, insert_size);
                }
            }
        }
        
        // Update alignment flags for paired-end information
        if let Some(ref mut a1) = alignment1 {
            self.update_paired_flags(a1, &alignment2, true, is_proper_pair);
        }
        
        if let Some(ref mut a2) = alignment2 {
            self.update_paired_flags(a2, &alignment1, false, is_proper_pair);
        }
        
        Ok(PairedAlignment {
            read1: alignment1,
            read2: alignment2,
            insert_size,
            is_proper_pair,
            rescued,
        })
    }
    
    /// Attempt mate rescue scanning for unmapped or poorly mapped reads
    fn attempt_mate_rescue(
        &self,
        aligner: &mut Aligner,
        read1: &Sequence,
        read2: &Sequence,
        alignment1: &Option<AlignmentResult>,
        alignment2: &Option<AlignmentResult>,
    ) -> Result<Option<(Option<AlignmentResult>, Option<AlignmentResult>)>> {
        
        // Determine which read is the anchor (well-mapped) and which needs rescue
        let (anchor_alignment, rescue_read, is_read1_anchor) = if let Some(ref a1) = alignment1 {
            if a1.mapq >= self.config.min_rescue_anchor_mapq {
                (a1, read2, true)
            } else if let Some(ref a2) = alignment2 {
                if a2.mapq >= self.config.min_rescue_anchor_mapq {
                    (a2, read1, false)
                } else {
                    return Ok(None); // No good anchor
                }
            } else {
                return Ok(None); // Only one alignment and it's poor quality
            }
        } else if let Some(ref a2) = alignment2 {
            if a2.mapq >= self.config.min_rescue_anchor_mapq {
                (a2, read1, false)
            } else {
                return Ok(None); // Only one alignment and it's poor quality
            }
        } else {
            return Ok(None); // No alignments to use as anchor
        };
        
        debug!("Attempting mate rescue: anchor at {}:{}, MAPQ={}", 
               anchor_alignment.reference_id, anchor_alignment.position, anchor_alignment.mapq);
        
        // Define rescue search region based on expected insert size
        let rescue_start = if anchor_alignment.position > self.config.max_rescue_distance {
            anchor_alignment.position - self.config.max_rescue_distance
        } else {
            0
        };
        let rescue_end = anchor_alignment.position + self.config.max_rescue_distance;
        
        // Perform targeted alignment in the rescue region
        if let Some(rescued_alignment) = self.rescue_align_in_region(
            aligner, 
            rescue_read, 
            &anchor_alignment.reference_id,
            rescue_start,
            rescue_end
        )? {
            debug!("Mate rescue successful: rescued read at {}:{}", 
                   rescued_alignment.reference_id, rescued_alignment.position);
            
            // Return updated alignments
            if is_read1_anchor {
                Ok(Some((alignment1.clone(), Some(rescued_alignment))))
            } else {
                Ok(Some((Some(rescued_alignment), alignment2.clone())))
            }
        } else {
            debug!("Mate rescue failed: no good alignment found in rescue region");
            Ok(None)
        }
    }
    
    /// Perform targeted alignment within a specific genomic region
    fn rescue_align_in_region(
        &self,
        aligner: &mut Aligner,
        read: &Sequence,
        reference_id: &str,
        start: u32,
        end: u32,
    ) -> Result<Option<AlignmentResult>> {
        debug!("Rescue aligning {} in region {}:{}-{}", read.id, reference_id, start, end);
        
        // Extract seeds from the read
        let seeds = aligner.extract_seeds(read)?;
        
        if seeds.is_empty() {
            return Ok(None);
        }
        
        // Find seed hits within the rescue region
        let mut rescue_hits = Vec::new();
        for (position, seed_seq, is_reverse) in seeds {
            if let Ok(hits) = aligner.find_seed_hits(&seed_seq, position, is_reverse) {
                for hit in hits {
                    // Filter hits to the rescue region
                    if hit.reference_position >= start && hit.reference_position <= end {
                        rescue_hits.push(hit);
                    }
                }
            }
        }
        
        debug!("Found {} seed hits in rescue region", rescue_hits.len());
        
        if rescue_hits.is_empty() {
            return Ok(None);
        }
        
        // Cluster and evaluate rescue hits
        let clustered_hits = aligner.cluster_hits(rescue_hits);
        let scored_candidates = aligner.evaluate_alignment_candidates(read, clustered_hits)?;
        
        if let Some((best_hit, best_score)) = scored_candidates.first() {
            // Extend the best rescue hit
            if let Some(mut alignment) = aligner.extend_alignment(read, best_hit)? {
                // Reduce MAPQ for rescued alignments to indicate lower confidence
                alignment.mapq = (alignment.mapq / 2).max(1);
                
                debug!("Rescue alignment: score={}, MAPQ={}", best_score, alignment.mapq);
                return Ok(Some(alignment));
            }
        }
        
        Ok(None)
    }
    
    /// Calculate insert size from two alignments
    fn calculate_insert_size(&self, alignment1: &AlignmentResult, alignment2: &AlignmentResult) -> Option<i32> {
        // Only calculate insert size for reads on the same reference
        if alignment1.reference_id != alignment2.reference_id {
            return None;
        }
        
        // Calculate insert size based on alignment positions and orientations
        let (left_pos, right_pos) = if alignment1.position <= alignment2.position {
            (alignment1.position, alignment2.position)
        } else {
            (alignment2.position, alignment1.position)
        };
        
        // Simple insert size calculation (can be refined)
        let insert_size = right_pos as i32 - left_pos as i32;
        
        // Adjust for read length to get template length
        let template_length = insert_size + 150; // Assuming ~150bp reads
        
        Some(template_length)
    }
    
    /// Check if two alignments form a proper pair
    fn is_proper_pair(&self, alignment1: &AlignmentResult, alignment2: &AlignmentResult, insert_size: Option<i32>) -> bool {
        // Must be on same reference
        if alignment1.reference_id != alignment2.reference_id {
            return false;
        }
        
        // Check insert size
        if let Some(isize) = insert_size {
            if !self.insert_stats.is_valid_insert_size(isize.abs()) {
                return false;
            }
        } else {
            return false;
        }
        
        // Check orientations (typical paired-end: R1 forward, R2 reverse)
        if alignment1.is_reverse == alignment2.is_reverse {
            return false; // Both same orientation
        }
        
        // Additional checks for proper pairing can be added here
        true
    }
    
    /// Update alignment flags for paired-end information
    fn update_paired_flags(
        &self,
        alignment: &mut AlignmentResult,
        mate_alignment: &Option<AlignmentResult>,
        is_first_in_pair: bool,
        is_proper_pair: bool,
    ) {
        alignment.is_paired = true;
        alignment.is_proper_pair = is_proper_pair;
        
        if let Some(mate) = mate_alignment {
            alignment.mate_reference_id = Some(mate.reference_id.clone());
            alignment.mate_position = Some(mate.position);
            
            // Calculate template length
            if alignment.reference_id == mate.reference_id {
                alignment.template_length = Some(mate.position as i32 - alignment.position as i32);
            }
        } else {
            alignment.mate_reference_id = Some("*".to_string());
            alignment.mate_position = Some(0);
            alignment.template_length = Some(0);
        }
    }
    
    /// Get current insert size statistics
    pub fn get_insert_stats(&self) -> &InsertSizeStats {
        &self.insert_stats
    }
    
    /// Update configuration with new insert size parameters
    pub fn update_insert_size_config(&mut self, mean: f64, std_dev: f64) {
        self.config.mean_insert_size = mean;
        self.config.insert_size_std = std_dev;
        debug!("Updated insert size config: mean={:.1}, std_dev={:.1}", mean, std_dev);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_insert_size_stats() {
        let mut stats = InsertSizeStats::new();
        
        // Add some observations
        for i in 400..600 {
            stats.observe(i);
        }
        
        stats.update_statistics();
        
        // Check that mean is reasonable
        assert!((stats.mean - 499.5).abs() < 1.0);
        assert!(stats.std_dev > 50.0 && stats.std_dev < 70.0);
        
        // Test validation
        assert!(stats.is_valid_insert_size(500));
        assert!(stats.is_valid_insert_size(300));
        assert!(stats.is_valid_insert_size(700));
        assert!(!stats.is_valid_insert_size(100));
        assert!(!stats.is_valid_insert_size(1000));
    }
    
    #[test]
    fn test_paired_end_config() {
        let config = PairedEndConfig::default();
        assert_eq!(config.min_insert_size, 50);
        assert_eq!(config.max_insert_size, 1000);
        assert_eq!(config.mean_insert_size, 500.0);
        assert!(config.enable_mate_rescue);
    }
    
    #[test]
    fn test_proper_pair_detection() {
        let config = PairedEndConfig::default();
        let mut pe_aligner = PairedEndAligner::new(config);
        
        let alignment1 = AlignmentResult {
            read_id: "read1".to_string(),
            reference_id: "chr1".to_string(),
            position: 1000,
            cigar: "100M".to_string(),
            mapq: 60,
            is_reverse: false,
            is_paired: false,
            is_proper_pair: false,
            mate_reference_id: None,
            mate_position: None,
            template_length: None,
        };
        
        let alignment2 = AlignmentResult {
            read_id: "read2".to_string(),
            reference_id: "chr1".to_string(),
            position: 1400,
            cigar: "100M".to_string(),
            mapq: 60,
            is_reverse: true,
            is_paired: false,
            is_proper_pair: false,
            mate_reference_id: None,
            mate_position: None,
            template_length: None,
        };
        
        let insert_size = pe_aligner.calculate_insert_size(&alignment1, &alignment2);
        assert!(insert_size.is_some());
        
        let is_proper = pe_aligner.is_proper_pair(&alignment1, &alignment2, insert_size);
        assert!(is_proper);
    }
}