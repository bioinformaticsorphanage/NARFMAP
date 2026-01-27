//! Smith-Waterman alignment using bio crate

use bio::alignment::pairwise::{self, Scoring};
use bio::alignment::AlignmentOperation;

use crate::alignment::mapper::{reverse_complement, CandidateRegion, SeedMapper};
use crate::reference::{HashTable, ReferenceSequence};
use crate::{Alignment, Read, ScoringParams};

/// Aligner for performing Smith-Waterman alignment
pub struct Aligner<'a> {
    hash_table: &'a HashTable,
    reference: &'a ReferenceSequence,
    scoring: ScoringParams,
    min_score: i32,
}

impl<'a> Aligner<'a> {
    /// Create a new aligner
    pub fn new(hash_table: &'a HashTable, reference: &'a ReferenceSequence) -> Self {
        Self {
            hash_table,
            reference,
            scoring: ScoringParams::default(),
            min_score: 22, // Default minimum alignment score
        }
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

    /// Align a read
    pub fn align(&self, read: &Read) -> Alignment {
        // Map read to find candidate regions
        let mapper = SeedMapper::new(self.hash_table);
        let candidates = mapper.map(read);

        if candidates.is_empty() {
            return Alignment::unmapped(read);
        }

        // Try to align to each candidate region
        let mut best_alignment: Option<Alignment> = None;
        let mut best_score = self.min_score - 1;

        for candidate in candidates.iter().take(5) {
            if let Some(aln) = self.align_to_region(read, candidate) {
                if aln.score > best_score {
                    best_score = aln.score;
                    best_alignment = Some(aln);
                }
            }
        }

        best_alignment.unwrap_or_else(|| Alignment::unmapped(read))
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
        let cigar = self.build_cigar(&alignment.operations);

        // Calculate mapping quality (simplified)
        let mapq = self.calculate_mapq(alignment.score, read.sequence.len());

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
            mapq,
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
    fn build_cigar(&self, ops: &[AlignmentOperation]) -> String {
        if ops.is_empty() {
            return "*".to_string();
        }

        let mut cigar = String::new();
        let mut current_op: Option<char> = None;
        let mut count = 0;

        for op in ops {
            let op_char = match op {
                AlignmentOperation::Match => 'M',
                AlignmentOperation::Subst => 'M', // Treat substitution as M (not X)
                AlignmentOperation::Ins => 'I',
                AlignmentOperation::Del => 'D',
                AlignmentOperation::Xclip(_) => 'S', // Soft clip
                AlignmentOperation::Yclip(_) => continue, // Skip reference clips
            };

            if current_op == Some(op_char) {
                count += 1;
            } else {
                if let Some(prev_op) = current_op {
                    cigar.push_str(&format!("{}{}", count, prev_op));
                }
                current_op = Some(op_char);
                count = 1;
            }
        }

        if let Some(op) = current_op {
            cigar.push_str(&format!("{}{}", count, op));
        }

        if cigar.is_empty() {
            "*".to_string()
        } else {
            cigar
        }
    }

    /// Calculate mapping quality
    fn calculate_mapq(&self, score: i32, read_len: usize) -> u8 {
        // Simplified MAPQ calculation
        // In reality, this should consider multiple alignments, etc.
        let max_score = (read_len as i32) * self.scoring.match_score;
        let score_ratio = (score as f64) / (max_score as f64);

        (score_ratio * 60.0).clamp(0.0, 60.0) as u8
    }
}
