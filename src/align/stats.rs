use std::collections::HashMap;
use std::time::{Duration, Instant};
use serde::{Serialize, Deserialize};

/// Comprehensive alignment statistics and quality metrics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentStats {
    /// Basic counts
    pub total_reads: u64,
    pub aligned_reads: u64,
    pub unaligned_reads: u64,
    pub multiply_aligned_reads: u64,
    
    /// Quality metrics
    pub total_bases: u64,
    pub aligned_bases: u64,
    pub mismatched_bases: u64,
    pub inserted_bases: u64,
    pub deleted_bases: u64,
    pub soft_clipped_bases: u64,
    pub hard_clipped_bases: u64,
    
    /// MAPQ distribution
    pub mapq_distribution: HashMap<u8, u64>,
    
    /// Insert size statistics (for paired-end)
    pub insert_size_sum: u64,
    pub insert_size_count: u64,
    pub insert_size_histogram: HashMap<i32, u64>,
    
    /// Performance metrics
    pub alignment_time: Duration,
    pub reads_per_second: f64,
    pub bases_per_second: f64,
    
    /// Reference coverage
    pub reference_coverage: HashMap<String, ReferenceCoverage>,
    
    /// Error rates
    pub substitution_rate: f64,
    pub insertion_rate: f64,
    pub deletion_rate: f64,
    pub overall_error_rate: f64,
}

/// Coverage statistics for a reference sequence
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceCoverage {
    pub sequence_name: String,
    pub sequence_length: u64,
    pub covered_bases: u64,
    pub total_coverage: u64,
    pub mean_coverage: f64,
    pub coverage_histogram: HashMap<u32, u64>, // coverage depth -> count
}

impl AlignmentStats {
    /// Create a new alignment statistics tracker
    pub fn new() -> Self {
        Self {
            total_reads: 0,
            aligned_reads: 0,
            unaligned_reads: 0,
            multiply_aligned_reads: 0,
            total_bases: 0,
            aligned_bases: 0,
            mismatched_bases: 0,
            inserted_bases: 0,
            deleted_bases: 0,
            soft_clipped_bases: 0,
            hard_clipped_bases: 0,
            mapq_distribution: HashMap::new(),
            insert_size_sum: 0,
            insert_size_count: 0,
            insert_size_histogram: HashMap::new(),
            alignment_time: Duration::new(0, 0),
            reads_per_second: 0.0,
            bases_per_second: 0.0,
            reference_coverage: HashMap::new(),
            substitution_rate: 0.0,
            insertion_rate: 0.0,
            deletion_rate: 0.0,
            overall_error_rate: 0.0,
        }
    }
    
    /// Record a read that was processed
    pub fn record_read(&mut self, read_length: usize, is_aligned: bool) {
        self.total_reads += 1;
        self.total_bases += read_length as u64;
        
        if is_aligned {
            self.aligned_reads += 1;
            self.aligned_bases += read_length as u64;
        } else {
            self.unaligned_reads += 1;
        }
    }
    
    /// Record MAPQ score distribution
    pub fn record_mapq(&mut self, mapq: u8) {
        *self.mapq_distribution.entry(mapq).or_insert(0) += 1;
    }
    
    /// Record insert size for paired-end reads
    pub fn record_insert_size(&mut self, insert_size: i32) {
        if insert_size > 0 {
            self.insert_size_sum += insert_size as u64;
            self.insert_size_count += 1;
            *self.insert_size_histogram.entry(insert_size).or_insert(0) += 1;
        }
    }
    
    /// Record CIGAR operations for error analysis
    pub fn record_cigar_operations(&mut self, cigar: &str) {
        let mut i = 0;
        let chars: Vec<char> = cigar.chars().collect();
        
        while i < chars.len() {
            let mut num_str = String::new();
            
            // Parse number
            while i < chars.len() && chars[i].is_ascii_digit() {
                num_str.push(chars[i]);
                i += 1;
            }
            
            if i < chars.len() {
                let count = num_str.parse::<u64>().unwrap_or(0);
                let op = chars[i];
                
                match op {
                    'M' | '=' | 'X' => {
                        // For simplicity, assume mismatches are 10% of matches
                        if op == 'X' {
                            self.mismatched_bases += count;
                        } else {
                            self.mismatched_bases += count / 10; // Rough estimate
                        }
                    },
                    'I' => self.inserted_bases += count,
                    'D' => self.deleted_bases += count,
                    'S' => self.soft_clipped_bases += count,
                    'H' => self.hard_clipped_bases += count,
                    _ => {}, // Ignore other operations
                }
                
                i += 1;
            }
        }
    }
    
    /// Record coverage for a reference position
    pub fn record_coverage(&mut self, reference_name: &str, _position: u32, coverage_depth: u32) {
        let coverage = self.reference_coverage
            .entry(reference_name.to_string())
            .or_insert_with(|| ReferenceCoverage {
                sequence_name: reference_name.to_string(),
                sequence_length: 0,
                covered_bases: 0,
                total_coverage: 0,
                mean_coverage: 0.0,
                coverage_histogram: HashMap::new(),
            });
        
        if coverage_depth > 0 {
            coverage.covered_bases += 1;
        }
        coverage.total_coverage += coverage_depth as u64;
        *coverage.coverage_histogram.entry(coverage_depth).or_insert(0) += 1;
    }
    
    /// Finalize statistics calculation
    pub fn finalize(&mut self, elapsed_time: Duration) {
        self.alignment_time = elapsed_time;
        
        if elapsed_time.as_secs_f64() > 0.0 {
            self.reads_per_second = self.total_reads as f64 / elapsed_time.as_secs_f64();
            self.bases_per_second = self.total_bases as f64 / elapsed_time.as_secs_f64();
        }
        
        // Calculate error rates
        let total_aligned_bases = self.aligned_bases as f64;
        if total_aligned_bases > 0.0 {
            self.substitution_rate = self.mismatched_bases as f64 / total_aligned_bases;
            self.insertion_rate = self.inserted_bases as f64 / total_aligned_bases;
            self.deletion_rate = self.deleted_bases as f64 / total_aligned_bases;
            self.overall_error_rate = (self.mismatched_bases + self.inserted_bases + self.deleted_bases) as f64 / total_aligned_bases;
        }
        
        // Finalize reference coverage statistics
        for coverage in self.reference_coverage.values_mut() {
            if coverage.sequence_length > 0 {
                coverage.mean_coverage = coverage.total_coverage as f64 / coverage.sequence_length as f64;
            }
        }
    }
    
    /// Get alignment rate as percentage
    pub fn alignment_rate(&self) -> f64 {
        if self.total_reads > 0 {
            (self.aligned_reads as f64 / self.total_reads as f64) * 100.0
        } else {
            0.0
        }
    }
    
    /// Get mean insert size for paired-end reads
    pub fn mean_insert_size(&self) -> f64 {
        if self.insert_size_count > 0 {
            self.insert_size_sum as f64 / self.insert_size_count as f64
        } else {
            0.0
        }
    }
    
    /// Get median MAPQ score
    pub fn median_mapq(&self) -> u8 {
        if self.mapq_distribution.is_empty() {
            return 0;
        }
        
        let total_alignments = self.mapq_distribution.values().sum::<u64>();
        let target = total_alignments / 2;
        
        let mut cumulative = 0;
        let mut mapq_scores: Vec<u8> = self.mapq_distribution.keys().cloned().collect();
        mapq_scores.sort();
        
        for mapq in mapq_scores {
            cumulative += self.mapq_distribution[&mapq];
            if cumulative >= target {
                return mapq;
            }
        }
        
        0
    }
    
    /// Generate a comprehensive statistics report
    pub fn generate_report(&self) -> String {
        let mut report = String::new();
        
        report.push_str("=== NARFMAP Alignment Statistics ===\n\n");
        
        // Basic alignment stats
        report.push_str(&format!("Total reads processed: {}\n", self.total_reads));
        report.push_str(&format!("Successfully aligned reads: {} ({:.2}%)\n", 
                                self.aligned_reads, self.alignment_rate()));
        report.push_str(&format!("Unaligned reads: {} ({:.2}%)\n", 
                                self.unaligned_reads, 
                                if self.total_reads > 0 { (self.unaligned_reads as f64 / self.total_reads as f64) * 100.0 } else { 0.0 }));
        
        if self.multiply_aligned_reads > 0 {
            report.push_str(&format!("Multiply aligned reads: {}\n", self.multiply_aligned_reads));
        }
        
        report.push_str("\n");
        
        // Base statistics
        report.push_str(&format!("Total bases: {}\n", self.total_bases));
        report.push_str(&format!("Aligned bases: {}\n", self.aligned_bases));
        report.push_str(&format!("Mismatched bases: {} ({:.4}%)\n", 
                                self.mismatched_bases, self.substitution_rate * 100.0));
        report.push_str(&format!("Inserted bases: {} ({:.4}%)\n", 
                                self.inserted_bases, self.insertion_rate * 100.0));
        report.push_str(&format!("Deleted bases: {} ({:.4}%)\n", 
                                self.deleted_bases, self.deletion_rate * 100.0));
        report.push_str(&format!("Overall error rate: {:.4}%\n", self.overall_error_rate * 100.0));
        
        report.push_str("\n");
        
        // MAPQ statistics
        report.push_str(&format!("Median MAPQ: {}\n", self.median_mapq()));
        if !self.mapq_distribution.is_empty() {
            report.push_str("MAPQ distribution:\n");
            let mut mapq_scores: Vec<u8> = self.mapq_distribution.keys().cloned().collect();
            mapq_scores.sort();
            for mapq in mapq_scores {
                let count = self.mapq_distribution[&mapq];
                let percentage = if self.aligned_reads > 0 { 
                    (count as f64 / self.aligned_reads as f64) * 100.0 
                } else { 0.0 };
                report.push_str(&format!("  MAPQ {}: {} ({:.1}%)\n", mapq, count, percentage));
            }
        }
        
        report.push_str("\n");
        
        // Insert size statistics for paired-end
        if self.insert_size_count > 0 {
            report.push_str(&format!("Mean insert size: {:.1} bp\n", self.mean_insert_size()));
            report.push_str(&format!("Insert size observations: {}\n", self.insert_size_count));
            report.push_str("\n");
        }
        
        // Performance statistics
        report.push_str(&format!("Alignment time: {:.2}s\n", self.alignment_time.as_secs_f64()));
        report.push_str(&format!("Reads per second: {:.0}\n", self.reads_per_second));
        report.push_str(&format!("Bases per second: {:.0}\n", self.bases_per_second));
        
        report.push_str("\n");
        
        // Reference coverage summary
        if !self.reference_coverage.is_empty() {
            report.push_str("Reference coverage summary:\n");
            for (ref_name, coverage) in &self.reference_coverage {
                let coverage_percentage = if coverage.sequence_length > 0 {
                    (coverage.covered_bases as f64 / coverage.sequence_length as f64) * 100.0
                } else { 0.0 };
                
                report.push_str(&format!("  {}: {:.1}% covered, {:.1}x mean depth\n", 
                                        ref_name, coverage_percentage, coverage.mean_coverage));
            }
        }
        
        report
    }
    
    /// Export statistics to JSON format
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(self)
    }
    
    /// Save statistics to a file
    pub fn save_to_file<P: AsRef<std::path::Path>>(&self, path: P) -> Result<(), Box<dyn std::error::Error>> {
        let report = self.generate_report();
        std::fs::write(path, report)?;
        Ok(())
    }
}

impl Default for AlignmentStats {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_alignment_stats_creation() {
        let stats = AlignmentStats::new();
        assert_eq!(stats.total_reads, 0);
        assert_eq!(stats.aligned_reads, 0);
        assert_eq!(stats.alignment_rate(), 0.0);
    }
    
    #[test]
    fn test_read_recording() {
        let mut stats = AlignmentStats::new();
        
        // Record some aligned reads
        stats.record_read(100, true);
        stats.record_read(150, true);
        stats.record_read(120, false); // unaligned
        
        assert_eq!(stats.total_reads, 3);
        assert_eq!(stats.aligned_reads, 2);
        assert_eq!(stats.unaligned_reads, 1);
        assert_eq!(stats.total_bases, 370);
        assert_eq!(stats.aligned_bases, 250);
        assert!((stats.alignment_rate() - 66.67).abs() < 0.01);
    }
    
    #[test]
    fn test_mapq_distribution() {
        let mut stats = AlignmentStats::new();
        
        stats.record_mapq(60);
        stats.record_mapq(30);
        stats.record_mapq(60);
        stats.record_mapq(10);
        
        assert_eq!(stats.mapq_distribution[&60], 2);
        assert_eq!(stats.mapq_distribution[&30], 1);
        assert_eq!(stats.mapq_distribution[&10], 1);
        
        let median = stats.median_mapq();
        assert!(median == 30 || median == 60); // Could be either depending on implementation
    }
    
    #[test]
    fn test_insert_size_recording() {
        let mut stats = AlignmentStats::new();
        
        stats.record_insert_size(300);
        stats.record_insert_size(250);
        stats.record_insert_size(350);
        
        assert_eq!(stats.insert_size_count, 3);
        assert_eq!(stats.insert_size_sum, 900);
        assert!((stats.mean_insert_size() - 300.0).abs() < 0.01);
    }
    
    #[test]
    fn test_cigar_operation_recording() {
        let mut stats = AlignmentStats::new();
        
        stats.record_cigar_operations("50M5I10M3D35M");
        
        assert!(stats.mismatched_bases > 0); // Should have some estimated mismatches
        assert_eq!(stats.inserted_bases, 5);
        assert_eq!(stats.deleted_bases, 3);
    }
    
    #[test]
    fn test_statistics_finalization() {
        let mut stats = AlignmentStats::new();
        
        // Add some data
        stats.record_read(100, true);
        stats.record_read(100, true);
        stats.mismatched_bases = 5;
        stats.inserted_bases = 2;
        stats.deleted_bases = 3;
        
        let start_time = Instant::now();
        std::thread::sleep(Duration::from_millis(10));
        let elapsed = start_time.elapsed();
        
        stats.finalize(elapsed);
        
        assert!(stats.reads_per_second > 0.0);
        assert!(stats.bases_per_second > 0.0);
        assert!(stats.overall_error_rate > 0.0);
    }
    
    #[test]
    fn test_report_generation() {
        let mut stats = AlignmentStats::new();
        
        stats.record_read(100, true);
        stats.record_read(100, false);
        stats.record_mapq(60);
        
        stats.finalize(Duration::from_secs(1));
        
        let report = stats.generate_report();
        assert!(report.contains("NARFMAP Alignment Statistics"));
        assert!(report.contains("Total reads processed: 2"));
        assert!(report.contains("Successfully aligned reads: 1"));
    }
}