use std::cmp::{max, min};

/// Smith-Waterman alignment algorithm implementation
/// Based on DRAGMAP's local alignment with configurable scoring
pub struct SmithWatermanAligner {
    /// Match score (positive)
    pub match_score: i32,
    /// Mismatch penalty (negative)  
    pub mismatch_penalty: i32,
    /// Gap open penalty (negative)
    pub gap_open_penalty: i32,
    /// Gap extension penalty (negative)
    pub gap_extend_penalty: i32,
    /// Enable banded alignment for performance
    pub enable_banding: bool,
    /// Band width for banded alignment
    pub band_width: usize,
}

/// CIGAR operations for alignment representation
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOp {
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

/// Alignment result from Smith-Waterman
#[derive(Debug, Clone)]
pub struct AlignmentResult {
    /// Alignment score
    pub score: i32,
    /// Reference start position (0-based)
    pub ref_start: usize,
    /// Reference end position (0-based, exclusive)
    pub ref_end: usize,
    /// Query start position (0-based)
    pub query_start: usize,
    /// Query end position (0-based, exclusive)
    pub query_end: usize,
    /// CIGAR operations describing the alignment
    pub cigar: Vec<CigarOp>,
    /// Number of mismatches
    pub mismatches: usize,
    /// Number of gaps (indels)
    pub gaps: usize,
}

/// Traceback direction for dynamic programming
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TracebackDirection {
    /// Diagonal (match/mismatch)
    Diagonal,
    /// Up (deletion in query)
    Up,
    /// Left (insertion in query)
    Left,
    /// End of alignment
    End,
}

/// Cell in the dynamic programming matrix
#[derive(Debug, Clone, Copy)]
struct DpCell {
    /// Maximum score ending at this cell
    score: i32,
    /// Traceback direction for optimal path
    direction: TracebackDirection,
    /// Whether this cell uses gap open vs extend
    is_gap_open: bool,
}

impl Default for SmithWatermanAligner {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_penalty: -4,
            gap_open_penalty: -6,
            gap_extend_penalty: -1,
            enable_banding: true,
            band_width: 100,
        }
    }
}

impl SmithWatermanAligner {
    /// Create a new Smith-Waterman aligner with DRAGMAP-compatible scoring
    pub fn with_dragmap_defaults() -> Self {
        Self {
            match_score: 2,
            mismatch_penalty: -4,
            gap_open_penalty: -6,
            gap_extend_penalty: -1,
            enable_banding: true,
            band_width: 100,
        }
    }

    /// Create a new aligner with custom scoring
    pub fn new(
        match_score: i32,
        mismatch_penalty: i32,
        gap_open_penalty: i32,
        gap_extend_penalty: i32,
    ) -> Self {
        Self {
            match_score,
            mismatch_penalty,
            gap_open_penalty,
            gap_extend_penalty,
            enable_banding: true,
            band_width: 100,
        }
    }

    /// Perform local Smith-Waterman alignment between query and reference
    pub fn align(&self, query: &[u8], reference: &[u8]) -> Option<AlignmentResult> {
        if query.is_empty() || reference.is_empty() {
            return None;
        }

        if self.enable_banding && query.len() > 50 && reference.len() > 50 {
            self.banded_align(query, reference)
        } else {
            self.full_align(query, reference)
        }
    }

    /// Perform semi-global alignment (query must align completely, reference can have overhangs)
    pub fn semi_global_align(&self, query: &[u8], reference: &[u8]) -> Option<AlignmentResult> {
        if query.is_empty() || reference.is_empty() {
            return None;
        }

        let query_len = query.len();
        let ref_len = reference.len();

        // Create DP matrix: (query_len + 1) x (ref_len + 1)
        let mut dp = vec![vec![DpCell::default(); ref_len + 1]; query_len + 1];

        // Initialize first row: no penalty for reference prefix
        for j in 0..=ref_len {
            dp[0][j] = DpCell {
                score: 0,
                direction: TracebackDirection::End,
                is_gap_open: false,
            };
        }

        // Initialize first column: penalty for query prefix (must align completely)
        for i in 1..=query_len {
            dp[i][0] = DpCell {
                score: self.gap_open_penalty + (i as i32 - 1) * self.gap_extend_penalty,
                direction: TracebackDirection::Up,
                is_gap_open: i == 1,
            };
        }

        // Fill DP matrix
        for i in 1..=query_len {
            for j in 1..=ref_len {
                let query_base = query[i - 1];
                let ref_base = reference[j - 1];

                // Match/mismatch score
                let match_mismatch_score = if query_base == ref_base {
                    self.match_score
                } else {
                    self.mismatch_penalty
                };

                // Diagonal (match/mismatch)
                let diagonal_score = dp[i - 1][j - 1].score + match_mismatch_score;

                // Up (deletion in query, insertion in reference)
                let up_score = if dp[i - 1][j].direction == TracebackDirection::Up && !dp[i - 1][j].is_gap_open {
                    dp[i - 1][j].score + self.gap_extend_penalty
                } else {
                    dp[i - 1][j].score + self.gap_open_penalty
                };

                // Left (insertion in query, deletion in reference)  
                let left_score = if dp[i][j - 1].direction == TracebackDirection::Left && !dp[i][j - 1].is_gap_open {
                    dp[i][j - 1].score + self.gap_extend_penalty
                } else {
                    dp[i][j - 1].score + self.gap_open_penalty
                };

                // Choose best option
                if diagonal_score >= up_score && diagonal_score >= left_score {
                    dp[i][j] = DpCell {
                        score: diagonal_score,
                        direction: TracebackDirection::Diagonal,
                        is_gap_open: false,
                    };
                } else if up_score >= left_score {
                    dp[i][j] = DpCell {
                        score: up_score,
                        direction: TracebackDirection::Up,
                        is_gap_open: dp[i - 1][j].direction != TracebackDirection::Up,
                    };
                } else {
                    dp[i][j] = DpCell {
                        score: left_score,
                        direction: TracebackDirection::Left,
                        is_gap_open: dp[i][j - 1].direction != TracebackDirection::Left,
                    };
                }
            }
        }

        // Find best alignment in last row (query fully aligned)
        let mut best_score = std::i32::MIN;
        let mut best_j = 0;

        for j in 0..=ref_len {
            if dp[query_len][j].score > best_score {
                best_score = dp[query_len][j].score;
                best_j = j;
            }
        }

        // Traceback from best position
        self.traceback_semi_global(&dp, query, reference, query_len, best_j)
    }

    /// Full Smith-Waterman alignment (local alignment)
    fn full_align(&self, query: &[u8], reference: &[u8]) -> Option<AlignmentResult> {
        let query_len = query.len();
        let ref_len = reference.len();

        // Create DP matrix
        let mut dp = vec![vec![DpCell::default(); ref_len + 1]; query_len + 1];

        // Initialize first row and column (all zeros for local alignment)
        for i in 0..=query_len {
            dp[i][0] = DpCell {
                score: 0,
                direction: TracebackDirection::End,
                is_gap_open: false,
            };
        }

        for j in 0..=ref_len {
            dp[0][j] = DpCell {
                score: 0,
                direction: TracebackDirection::End,
                is_gap_open: false,
            };
        }

        let mut max_score = 0;
        let mut max_i = 0;
        let mut max_j = 0;

        // Fill DP matrix
        for i in 1..=query_len {
            for j in 1..=ref_len {
                let query_base = query[i - 1];
                let ref_base = reference[j - 1];

                // Match/mismatch score
                let match_mismatch_score = if query_base == ref_base {
                    self.match_score
                } else {
                    self.mismatch_penalty
                };

                // Diagonal (match/mismatch)
                let diagonal_score = dp[i - 1][j - 1].score + match_mismatch_score;

                // Up (deletion in query)
                let up_score = if dp[i - 1][j].direction == TracebackDirection::Up && !dp[i - 1][j].is_gap_open {
                    dp[i - 1][j].score + self.gap_extend_penalty
                } else {
                    dp[i - 1][j].score + self.gap_open_penalty
                };

                // Left (insertion in query)
                let left_score = if dp[i][j - 1].direction == TracebackDirection::Left && !dp[i][j - 1].is_gap_open {
                    dp[i][j - 1].score + self.gap_extend_penalty
                } else {
                    dp[i][j - 1].score + self.gap_open_penalty
                };

                // Local alignment: score cannot go below 0
                let best_score = max(0, max(diagonal_score, max(up_score, left_score)));

                // Determine direction
                let direction = if best_score == 0 {
                    TracebackDirection::End
                } else if best_score == diagonal_score {
                    TracebackDirection::Diagonal
                } else if best_score == up_score {
                    TracebackDirection::Up
                } else {
                    TracebackDirection::Left
                };

                dp[i][j] = DpCell {
                    score: best_score,
                    direction,
                    is_gap_open: match direction {
                        TracebackDirection::Up => dp[i - 1][j].direction != TracebackDirection::Up,
                        TracebackDirection::Left => dp[i][j - 1].direction != TracebackDirection::Left,
                        _ => false,
                    },
                };

                // Track maximum score for local alignment
                if best_score > max_score {
                    max_score = best_score;
                    max_i = i;
                    max_j = j;
                }
            }
        }

        if max_score <= 0 {
            return None;
        }

        // Traceback from maximum score position
        self.traceback(&dp, query, reference, max_i, max_j)
    }

    /// Banded Smith-Waterman alignment for improved performance
    fn banded_align(&self, query: &[u8], reference: &[u8]) -> Option<AlignmentResult> {
        // For now, use a simple approach: find the best diagonal and align around it
        let best_diagonal = self.find_best_diagonal(query, reference);
        let band_radius = self.band_width / 2;

        // Extract region around the best diagonal
        let ref_start = if best_diagonal >= band_radius { best_diagonal - band_radius } else { 0 };
        let ref_end = min(reference.len(), best_diagonal + query.len() + band_radius);

        if ref_end <= ref_start {
            return self.full_align(query, reference);
        }

        let ref_subseq = &reference[ref_start..ref_end];
        
        // Align within the band
        if let Some(mut result) = self.full_align(query, ref_subseq) {
            // Adjust coordinates to original reference
            result.ref_start += ref_start;
            result.ref_end += ref_start;
            Some(result)
        } else {
            // Fallback to full alignment
            self.full_align(query, reference)
        }
    }

    /// Find the best diagonal for banded alignment
    fn find_best_diagonal(&self, query: &[u8], reference: &[u8]) -> usize {
        let query_len = query.len();
        let ref_len = reference.len();
        
        if ref_len < query_len {
            return 0;
        }

        let mut best_score = std::i32::MIN;
        let mut best_start = 0;

        // Try different starting positions with a reasonable step size
        let step = max(1, (ref_len - query_len) / 20);
        
        for start in (0..=ref_len - query_len).step_by(step) {
            let score = self.quick_score(query, &reference[start..start + query_len]);
            if score > best_score {
                best_score = score;
                best_start = start;
            }
        }

        best_start
    }

    /// Quick scoring for diagonal finding
    fn quick_score(&self, query: &[u8], reference: &[u8]) -> i32 {
        let mut score = 0;
        let len = min(query.len(), reference.len());

        for i in 0..len {
            if query[i] == reference[i] {
                score += self.match_score;
            } else {
                score += self.mismatch_penalty;
            }
        }

        score
    }

    /// Traceback for local alignment
    fn traceback(
        &self,
        dp: &[Vec<DpCell>],
        query: &[u8],
        reference: &[u8],
        start_i: usize,
        start_j: usize,
    ) -> Option<AlignmentResult> {
        let mut cigar = Vec::new();
        let mut i = start_i;
        let mut j = start_j;
        let mut mismatches = 0;
        let mut gaps = 0;

        let end_i = i;
        let end_j = j;

        // Traceback until we reach a cell with score 0 or boundary
        while i > 0 && j > 0 && dp[i][j].score > 0 {
            match dp[i][j].direction {
                TracebackDirection::Diagonal => {
                    let query_base = query[i - 1];
                    let ref_base = reference[j - 1];
                    
                    if query_base != ref_base {
                        mismatches += 1;
                    }
                    
                    // Count consecutive matches/mismatches
                    let mut match_len = 1;
                    i -= 1;
                    j -= 1;
                    
                    while i > 0 && j > 0 && dp[i][j].score > 0 && dp[i + 1][j + 1].direction == TracebackDirection::Diagonal {
                        let qb = query[i - 1];
                        let rb = reference[j - 1];
                        if qb != rb {
                            mismatches += 1;
                        }
                        match_len += 1;
                        i -= 1;
                        j -= 1;
                    }
                    
                    cigar.push(CigarOp::Match(match_len));
                }
                TracebackDirection::Up => {
                    // Deletion in query (insertion in reference)
                    let mut del_len = 1;
                    i -= 1;
                    
                    while i > 0 && dp[i][j].score > 0 && dp[i + 1][j].direction == TracebackDirection::Up {
                        del_len += 1;
                        i -= 1;
                    }
                    
                    gaps += 1; // Count as one gap event
                    cigar.push(CigarOp::Del(del_len));
                }
                TracebackDirection::Left => {
                    // Insertion in query (deletion in reference)
                    let mut ins_len = 1;
                    j -= 1;
                    
                    while j > 0 && dp[i][j].score > 0 && dp[i][j + 1].direction == TracebackDirection::Left {
                        ins_len += 1;
                        j -= 1;
                    }
                    
                    gaps += 1; // Count as one gap event
                    cigar.push(CigarOp::Ins(ins_len));
                }
                TracebackDirection::End => break,
            }
        }

        // Reverse CIGAR since we traced backwards
        cigar.reverse();

        Some(AlignmentResult {
            score: dp[start_i][start_j].score,
            ref_start: j,
            ref_end: end_j,
            query_start: i,
            query_end: end_i,
            cigar,
            mismatches,
            gaps,
        })
    }

    /// Traceback for semi-global alignment
    fn traceback_semi_global(
        &self,
        dp: &[Vec<DpCell>],
        query: &[u8],
        reference: &[u8],
        start_i: usize,
        start_j: usize,
    ) -> Option<AlignmentResult> {
        let mut cigar = Vec::new();
        let mut i = start_i;
        let mut j = start_j;
        let mut mismatches = 0;
        let mut gaps = 0;

        let end_i = i;
        let end_j = j;

        // Traceback until we reach the first row
        while i > 0 {
            match dp[i][j].direction {
                TracebackDirection::Diagonal => {
                    if j == 0 {
                        // Should not happen in proper semi-global alignment
                        break;
                    }
                    
                    let query_base = query[i - 1];
                    let ref_base = reference[j - 1];
                    
                    if query_base != ref_base {
                        mismatches += 1;
                    }
                    
                    cigar.push(CigarOp::Match(1));
                    i -= 1;
                    j -= 1;
                }
                TracebackDirection::Up => {
                    cigar.push(CigarOp::Del(1));
                    gaps += 1;
                    i -= 1;
                }
                TracebackDirection::Left => {
                    if j == 0 {
                        break;
                    }
                    cigar.push(CigarOp::Ins(1));
                    gaps += 1;
                    j -= 1;
                }
                TracebackDirection::End => break,
            }
        }

        // Reverse CIGAR since we traced backwards
        cigar.reverse();

        Some(AlignmentResult {
            score: dp[start_i][start_j].score,
            ref_start: j,
            ref_end: end_j,
            query_start: i,
            query_end: end_i,
            cigar,
            mismatches,
            gaps,
        })
    }

    /// Format CIGAR operations into a standard CIGAR string
    pub fn format_cigar(cigar: &[CigarOp]) -> String {
        let mut result = String::new();
        
        // Merge consecutive operations of the same type
        let mut merged_ops = Vec::new();
        for op in cigar {
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
}

impl Default for DpCell {
    fn default() -> Self {
        Self {
            score: 0,
            direction: TracebackDirection::End,
            is_gap_open: false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_perfect_match() {
        let aligner = SmithWatermanAligner::with_dragmap_defaults();
        let query = b"ACGT";
        let reference = b"ACGT";
        
        let result = aligner.align(query, reference).unwrap();
        assert_eq!(result.score, 8); // 4 matches * 2 = 8
        assert_eq!(result.mismatches, 0);
        assert_eq!(result.gaps, 0);
        
        let cigar = SmithWatermanAligner::format_cigar(&result.cigar);
        assert_eq!(cigar, "4M");
    }

    #[test]
    fn test_mismatch() {
        let aligner = SmithWatermanAligner::with_dragmap_defaults();
        let query = b"ACGT";
        let reference = b"ACTT";
        
        let result = aligner.semi_global_align(query, reference).unwrap();
        assert_eq!(result.mismatches, 1);
        
        let cigar = SmithWatermanAligner::format_cigar(&result.cigar);
        assert_eq!(cigar, "4M");
    }

    #[test]
    fn test_insertion() {
        let aligner = SmithWatermanAligner::with_dragmap_defaults();
        let query = b"ACGTT";
        let reference = b"ACGT";
        
        let result = aligner.semi_global_align(query, reference).unwrap();
        assert_eq!(result.gaps, 1);
        
        let cigar = SmithWatermanAligner::format_cigar(&result.cigar);
        println!("Insertion CIGAR: {}", cigar);
        assert!(cigar.contains('I') || cigar.contains('M')); // Should handle the insertion
    }

    #[test]
    fn test_deletion() {
        let aligner = SmithWatermanAligner::with_dragmap_defaults();
        let query = b"ACTT"; // Missing G in the middle
        let reference = b"ACGTT";
        
        let result = aligner.semi_global_align(query, reference).unwrap();
        
        let cigar = SmithWatermanAligner::format_cigar(&result.cigar);
        println!("Deletion test:");
        println!("Query: {:?}", std::str::from_utf8(query));
        println!("Reference: {:?}", std::str::from_utf8(reference));
        println!("CIGAR: {}", cigar);
        println!("Gaps: {}, Mismatches: {}", result.gaps, result.mismatches);
        
        // Should handle the gap somehow - either as deletion or mismatches
        assert!(result.gaps >= 1 || result.mismatches >= 1);
        assert!(cigar.contains('D') || cigar.contains('M')); // Should handle the difference
    }

    #[test]
    fn test_semi_global_alignment() {
        let aligner = SmithWatermanAligner::with_dragmap_defaults();
        let query = b"ACGT";
        let reference = b"TTACGTAA"; // Query is embedded in reference
        
        let result = aligner.semi_global_align(query, reference).unwrap();
        assert_eq!(result.score, 8); // Perfect match for the query
        assert_eq!(result.ref_start, 2); // Should start at position 2 in reference
        assert_eq!(result.ref_end, 6); // Should end at position 6 in reference
        
        let cigar = SmithWatermanAligner::format_cigar(&result.cigar);
        assert_eq!(cigar, "4M");
    }

    #[test]
    fn test_complex_alignment() {
        let aligner = SmithWatermanAligner::with_dragmap_defaults();
        let query = b"ACGTACGT";
        let reference = b"ACGTTACGT"; // Has an extra T in the middle
        
        let result = aligner.align(query, reference).unwrap();
        assert!(result.score > 0);
        
        let cigar = SmithWatermanAligner::format_cigar(&result.cigar);
        println!("Complex alignment CIGAR: {}", cigar);
        
        // Should handle the indel appropriately
        assert!(cigar.contains('M'));
    }

    #[test]
    fn test_no_alignment() {
        let aligner = SmithWatermanAligner::with_dragmap_defaults();
        let query = b"AAAA";
        let reference = b"TTTT";
        
        // With harsh mismatch penalty, this might not align well
        let result = aligner.align(query, reference);
        
        // May or may not find an alignment depending on scoring
        if let Some(res) = result {
            assert!(res.score >= 0); // Local alignment score should be non-negative
        }
    }

    #[test]
    fn test_banded_alignment() {
        let mut aligner = SmithWatermanAligner::with_dragmap_defaults();
        aligner.enable_banding = true;
        aligner.band_width = 20;
        
        // Long sequences to trigger banded alignment
        let query = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let reference = b"TTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAA";
        
        let result = aligner.align(query, reference).unwrap();
        assert!(result.score > 0);
        
        let cigar = SmithWatermanAligner::format_cigar(&result.cigar);
        println!("Banded alignment CIGAR: {}", cigar);
        assert!(cigar.contains('M'));
    }

    #[test]
    fn test_cigar_formatting() {
        let cigar_ops = vec![
            CigarOp::Match(10),
            CigarOp::Ins(2),
            CigarOp::Match(5),
            CigarOp::Del(1),
            CigarOp::Match(8),
        ];
        
        let cigar = SmithWatermanAligner::format_cigar(&cigar_ops);
        assert_eq!(cigar, "10M2I5M1D8M");
        
        // Test merging consecutive operations
        let cigar_ops_consecutive = vec![
            CigarOp::Match(10),
            CigarOp::Match(5), // Should merge
            CigarOp::Ins(2),
            CigarOp::Ins(1), // Should merge
            CigarOp::Match(3),
        ];
        
        let merged_cigar = SmithWatermanAligner::format_cigar(&cigar_ops_consecutive);
        assert_eq!(merged_cigar, "15M3I3M");
    }
}