//! Smith-Waterman alignment using bio crate

use bio::alignment::pairwise::{self, Scoring};
use bio::alignment::AlignmentOperation;

use crate::alignment::mapper::{reverse_complement, CandidateRegion, SeedMapper};
use crate::alignment::mapq::{compute_mapq, MAPQ_MAX};
use crate::reference::{HashTable, ReferenceSequence};
use crate::{Alignment, Read, ScoringParams};

/// Configuration for secondary alignment output
#[derive(Debug, Clone, Default)]
pub struct SecondaryAlignmentConfig {
    /// Maximum number of secondary alignments to output
    pub max_secondary: usize,
    /// Minimum score delta from primary to report secondary
    pub min_score_delta: i32,
    /// Minimum PHRED score delta for secondary
    pub min_phred_delta: i32,
}

/// Aligner for performing Smith-Waterman alignment
pub struct Aligner<'a> {
    hash_table: &'a HashTable,
    reference: &'a ReferenceSequence,
    scoring: ScoringParams,
    min_score: i32,
    secondary_config: SecondaryAlignmentConfig,
}

impl<'a> Aligner<'a> {
    /// Create a new aligner
    pub fn new(hash_table: &'a HashTable, reference: &'a ReferenceSequence) -> Self {
        Self {
            hash_table,
            reference,
            scoring: ScoringParams::default(),
            min_score: 22, // Default minimum alignment score
            secondary_config: SecondaryAlignmentConfig::default(),
        }
    }

    /// Configure secondary alignment output
    pub fn with_secondary_config(mut self, config: SecondaryAlignmentConfig) -> Self {
        self.secondary_config = config;
        self
    }

    /// Set scoring parameters
    pub fn with_scoring(mut self, scoring: ScoringParams) -> Self {
        self.scoring = scoring;
        self
    }

    /// Set minimum alignment score
    pub fn with_min_score(mut self, min_score: i32) -> Self {
        self.min_score = min_score;
        self
    }

    /// Align a read (returns primary alignment only)
    pub fn align(&self, read: &Read) -> Alignment {
        let alignments = self.align_all(read);
        alignments
            .into_iter()
            .next()
            .unwrap_or_else(|| Alignment::unmapped(read))
    }

    /// Align a read and return all alignments (primary + secondary)
    pub fn align_all(&self, read: &Read) -> Vec<Alignment> {
        // Map read to find candidate regions
        let mapper = SeedMapper::new(self.hash_table);
        let candidates = mapper.map(read);

        if candidates.is_empty() {
            return vec![Alignment::unmapped(read)];
        }

        // Try to align to each candidate region, collect all valid alignments
        let mut alignments: Vec<Alignment> = Vec::new();

        for candidate in candidates.iter().take(10) {
            if let Some(aln) = self.align_to_region(read, candidate) {
                alignments.push(aln);
            }
        }

        if alignments.is_empty() {
            return vec![Alignment::unmapped(read)];
        }

        // Sort by score descending
        alignments.sort_by(|a, b| b.score.cmp(&a.score));

        // Calculate proper MAPQ using best vs second-best scores
        let best_score = alignments[0].score;
        let second_score = alignments.get(1).map_or(0, |a| a.score);
        let snp_cost = self.scoring.mismatch_score.unsigned_abs() as i32;
        let mapq = compute_mapq(
            snp_cost.max(1),
            best_score,
            second_score,
            read.sequence.len() as i32,
        );

        // Update primary alignment MAPQ
        alignments[0].mapq = mapq.clamp(0, MAPQ_MAX) as u8;

        // Mark and filter secondary alignments
        let max_secondary = self.secondary_config.max_secondary;
        if max_secondary > 0 && alignments.len() > 1 {
            for aln in alignments.iter_mut().skip(1).take(max_secondary) {
                aln.flag |= 0x100; // Secondary alignment flag
                                   // Secondary alignments get MAPQ 0
                aln.mapq = 0;
            }
            // Keep primary + max_secondary
            alignments.truncate(1 + max_secondary);
        } else {
            // Only keep primary
            alignments.truncate(1);
        }

        alignments
    }

    /// Align read to a specific candidate region
    fn align_to_region(&self, read: &Read, region: &CandidateRegion) -> Option<Alignment> {
        // Get reference sequence for this region with some padding
        let padding = read.sequence.len() as u64;
        let ref_start = region.ref_start.saturating_sub(padding);
        let ref_len = (region.ref_end - region.ref_start + 2 * padding) as usize;
        let ref_seq = self.reference.get_sequence(ref_start, ref_len);

        // Get query sequence (reverse complement if needed)
        let query = if region.is_reverse {
            reverse_complement(&read.sequence)
        } else {
            read.sequence.clone()
        };

        // Create scoring matrix
        let scoring = Scoring::from_scores(
            self.scoring.gap_open,
            self.scoring.gap_extend,
            self.scoring.match_score,
            self.scoring.mismatch_score,
        )
        .xclip(0)
        .yclip(0);

        // Perform local alignment
        let mut aligner = pairwise::Aligner::with_scoring(scoring);
        let alignment = aligner.local(&query, &ref_seq);

        if alignment.score < self.min_score {
            return None;
        }

        // Build CIGAR string
        let cigar = Self::build_cigar(&alignment.operations);

        // Calculate alignment position
        let align_pos = ref_start + alignment.ystart as u64 + 1; // 1-based

        // Determine flags
        let mut flag = 0u16;
        if region.is_reverse {
            flag |= 0x10; // SEQ is reverse complemented
        }

        // Get reference name
        let ref_name = self
            .reference
            .find_sequence(align_pos)
            .map(|(seq, _)| seq.name.clone())
            .unwrap_or_else(|| "ref".to_string());

        Some(Alignment {
            read_name: read.name.clone(),
            flag,
            ref_name,
            position: align_pos,
            mapq: 0, // Will be calculated in align_all() based on best/second-best scores
            cigar,
            mate_ref_name: "*".to_string(),
            mate_position: 0,
            template_length: 0,
            sequence: read.sequence.clone(),
            quality: read.quality.clone(),
            score: alignment.score,
        })
    }

    /// Build CIGAR string from alignment operations
    fn build_cigar(ops: &[AlignmentOperation]) -> String {
        use std::fmt::Write;

        if ops.is_empty() {
            return "*".to_string();
        }

        let mut cigar = String::new();
        let mut current_op: Option<char> = None;
        let mut count = 0;

        for op in ops {
            let op_char = match op {
                AlignmentOperation::Match | AlignmentOperation::Subst => 'M',
                AlignmentOperation::Ins => 'I',
                AlignmentOperation::Del => 'D',
                AlignmentOperation::Xclip(_) => 'S', // Soft clip
                AlignmentOperation::Yclip(_) => continue, // Skip reference clips
            };

            if current_op == Some(op_char) {
                count += 1;
            } else {
                if let Some(prev_op) = current_op {
                    let _ = write!(cigar, "{count}{prev_op}");
                }
                current_op = Some(op_char);
                count = 1;
            }
        }

        if let Some(op) = current_op {
            let _ = write!(cigar, "{count}{op}");
        }

        if cigar.is_empty() {
            "*".to_string()
        } else {
            cigar
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::hash_table::HashTable;
    use crate::reference::ReferenceSequence;
    use std::path::PathBuf;

    fn tiny_ref_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("data")
            .join("tiny")
            .join("tiny-2x1Xrepeats.v8")
    }

    #[test]
    fn build_cigar_collapses_operations() {
        let ops = vec![
            AlignmentOperation::Match,
            AlignmentOperation::Match,
            AlignmentOperation::Ins,
            AlignmentOperation::Ins,
            AlignmentOperation::Del,
            AlignmentOperation::Subst,
            AlignmentOperation::Subst,
        ];

        let cigar = Aligner::build_cigar(&ops);

        assert_eq!(cigar, "2M2I1D2M");
    }

    #[test]
    fn build_cigar_skips_yclip_only() {
        let ops = vec![AlignmentOperation::Yclip(3)];

        let cigar = Aligner::build_cigar(&ops);

        assert_eq!(cigar, "*");
    }

    #[test]
    fn align_all_returns_primary_with_mapq() {
        let ref_dir = tiny_ref_dir();
        assert!(ref_dir.exists(), "missing tiny reference data");

        let hash_table = HashTable::load(&ref_dir).unwrap();
        let reference = ReferenceSequence::load(&ref_dir).unwrap();
        let aligner = Aligner::new(&hash_table, &reference);

        let read = Read {
            name: "test".to_string(),
            sequence: reference.get_sequence(163_840 + 71, 50),
            quality: vec![b'I'; 50],
        };

        let alignments = aligner.align_all(&read);
        assert!(!alignments.is_empty());
        // Primary alignment should have MAPQ set
        assert!(alignments[0].mapq > 0 || alignments[0].flag & 4 != 0);
    }
}
