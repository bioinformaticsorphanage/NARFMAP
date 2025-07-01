use std::path::{Path, PathBuf};
use tempfile::TempDir;
use assert_fs::prelude::*;

/// Test data paths
pub fn get_test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data")
}

pub fn get_tiny_fasta() -> PathBuf {
    get_test_data_dir().join("tiny/tiny-2x1Xrepeats.v8/tiny.fasta")
}

pub fn get_tiny_fastq() -> PathBuf {
    get_test_data_dir().join("tiny/1read.fastq")
}

pub fn get_paired_fastq() -> (PathBuf, PathBuf) {
    let data_dir = get_test_data_dir().join("tiny");
    (
        data_dir.join("1read-1X.fastq"),
        data_dir.join("1read-1X-r2.fastq")
    )
}

/// Create a temporary workspace for testing
pub fn create_test_workspace() -> TempDir {
    tempfile::tempdir().expect("Failed to create temporary directory")
}

/// Helper to verify SAM output format
pub fn verify_sam_format(sam_file: &Path) -> Result<(), Box<dyn std::error::Error>> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    
    let file = File::open(sam_file)?;
    let reader = BufReader::new(file);
    
    let mut has_header = false;
    let mut has_alignment = false;
    
    for line_result in reader.lines() {
        let line = line_result?;
        if line.starts_with("@") {
            has_header = true;
        } else if !line.trim().is_empty() {
            // Basic SAM format validation (11+ columns)
            let fields: Vec<&str> = line.split('\t').collect();
            assert!(fields.len() >= 11, "SAM line should have at least 11 fields: {}", line);
            has_alignment = true;
        }
    }
    
    assert!(has_header, "SAM file should contain header lines");
    Ok(())
}

/// Helper to count reads/alignments in SAM file
pub fn count_sam_records(sam_file: &Path) -> Result<(usize, usize), Box<dyn std::error::Error>> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    
    let file = File::open(sam_file)?;
    let reader = BufReader::new(file);
    
    let mut header_count = 0;
    let mut alignment_count = 0;
    
    for line_result in reader.lines() {
        let line = line_result?;
        if line.starts_with("@") {
            header_count += 1;
        } else if !line.trim().is_empty() {
            alignment_count += 1;
        }
    }
    
    Ok((header_count, alignment_count))
}

/// Verify test data exists
pub fn check_test_data_exists() -> Result<(), Box<dyn std::error::Error>> {
    let tiny_fasta = get_tiny_fasta();
    let tiny_fastq = get_tiny_fastq();
    let (fastq1, fastq2) = get_paired_fastq();
    
    if !tiny_fasta.exists() {
        return Err(format!("Test FASTA file not found: {}", tiny_fasta.display()).into());
    }
    
    if !tiny_fastq.exists() {
        return Err(format!("Test FASTQ file not found: {}", tiny_fastq.display()).into());
    }
    
    if !fastq1.exists() {
        return Err(format!("Test paired FASTQ file 1 not found: {}", fastq1.display()).into());
    }
    
    if !fastq2.exists() {
        return Err(format!("Test paired FASTQ file 2 not found: {}", fastq2.display()).into());
    }
    
    Ok(())
}

/// Helper to count lines in a file
pub fn count_file_lines(file_path: &Path) -> Result<usize, Box<dyn std::error::Error>> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let count = reader.lines().count();
    Ok(count)
}

/// Check if a file contains a specific string
pub fn file_contains_string(file_path: &Path, search_string: &str) -> Result<bool, Box<dyn std::error::Error>> {
    use std::fs;
    let contents = fs::read_to_string(file_path)?;
    Ok(contents.contains(search_string))
}