use assert_cmd::Command;
use assert_fs::prelude::*;
use predicates::prelude::*;
use test_case::test_case;

mod common;

#[test]
fn test_build_hashtable_basic() -> Result<(), Box<dyn std::error::Error>> {
    // Skip if test data doesn't exist
    if let Err(e) = common::check_test_data_exists() {
        println!("Skipping basic hashtable test due to missing test data: {}", e);
        return Ok(());
    }
    
    let temp_workspace = common::create_test_workspace();
    
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("build-hash-table")
        .arg("-r").arg(common::get_tiny_fasta())
        .arg("-k").arg("21")
        .arg("-t").arg("2")
        .current_dir(&temp_workspace);
    
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("Hash table generation complete"));
    
    // Check that hash table directory was created
    let hash_dir = temp_workspace.child("hash_table");
    hash_dir.assert(predicate::path::exists());
    
    // Check for specific files
    let hash_bin = hash_dir.child("hash_table.bin");
    let hash_cfg = hash_dir.child("hash_table.cfg");
    
    // These files should exist after successful hash table generation
    if hash_bin.path().exists() {
        println!("Hash table binary size: {} bytes", 
                 std::fs::metadata(hash_bin.path())?.len());
    }
    
    if hash_cfg.path().exists() {
        println!("Hash table config size: {} bytes", 
                 std::fs::metadata(hash_cfg.path())?.len());
    }
    
    Ok(())
}

#[test]
fn test_build_hashtable_missing_reference() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("build-hash-table")
        .arg("-r").arg("nonexistent.fasta");
    
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("Failed to build hash table").or(
            predicate::str::contains("No such file")
        ));
    
    Ok(())
}

#[test_case(15)]
#[test_case(19)] 
#[test_case(21)]
#[test_case(25)]
fn test_build_hashtable_different_kmer_sizes(k: usize) -> Result<(), Box<dyn std::error::Error>> {
    // Skip if test data doesn't exist
    if let Err(e) = common::check_test_data_exists() {
        println!("Skipping k-mer size test for k={} due to missing test data: {}", k, e);
        return Ok(());
    }
    
    let temp_workspace = common::create_test_workspace();
    
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("build-hash-table")
        .arg("-r").arg(common::get_tiny_fasta())
        .arg("-k").arg(k.to_string())
        .arg("-t").arg("1")
        .current_dir(&temp_workspace);
    
    let result = cmd.assert().success();
    result.stdout(predicate::str::contains(format!("k-mer size: {}", k)));
    
    // Verify hash table was created
    let hash_dir = temp_workspace.child("hash_table");
    hash_dir.assert(predicate::path::exists());
    
    println!("Successfully built hash table with k={}", k);
    
    Ok(())
}

#[test]
fn test_build_hashtable_different_thread_counts() -> Result<(), Box<dyn std::error::Error>> {
    // Skip if test data doesn't exist
    if let Err(e) = common::check_test_data_exists() {
        println!("Skipping thread count test due to missing test data: {}", e);
        return Ok(());
    }
    
    for threads in [1, 2, 4] {
        let temp_workspace = common::create_test_workspace();
        
        let mut cmd = Command::cargo_bin("narfmap")?;
        cmd.arg("build-hash-table")
            .arg("-r").arg(common::get_tiny_fasta())
            .arg("-k").arg("15")
            .arg("-t").arg(threads.to_string())
            .current_dir(&temp_workspace);
        
        cmd.assert()
            .success()
            .stdout(predicate::str::contains(format!("Using {} threads", threads)));
        
        println!("Successfully built hash table with {} threads", threads);
    }
    
    Ok(())
}

#[test]
fn test_build_hashtable_output_verification() -> Result<(), Box<dyn std::error::Error>> {
    // Skip if test data doesn't exist
    if let Err(e) = common::check_test_data_exists() {
        println!("Skipping output verification test due to missing test data: {}", e);
        return Ok(());
    }
    
    let temp_workspace = common::create_test_workspace();
    
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("build-hash-table")
        .arg("-r").arg(common::get_tiny_fasta())
        .arg("-k").arg("15")
        .arg("-t").arg("1")
        .current_dir(&temp_workspace);
    
    let output = cmd.assert().success();
    
    // Check that the command provides useful output
    output.stdout(predicate::str::contains("Building hash table from reference"));
    output.stdout(predicate::str::contains("Output directory"));
    output.stdout(predicate::str::contains("Hash table generation complete"));
    
    // Verify the actual output files
    let hash_dir = temp_workspace.child("hash_table");
    hash_dir.assert(predicate::path::exists());
    
    // Count files in hash table directory
    let hash_dir_entries: Vec<_> = std::fs::read_dir(hash_dir.path())?
        .filter_map(|entry| entry.ok())
        .collect();
    
    println!("Hash table directory contains {} files", hash_dir_entries.len());
    for entry in &hash_dir_entries {
        println!("  - {}", entry.file_name().to_string_lossy());
    }
    
    // Should have at least some files
    assert!(hash_dir_entries.len() > 0, "Hash table directory should contain files");
    
    Ok(())
}

#[test]
fn test_build_hashtable_decompression_error() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("build-hash-table")
        .arg("-r").arg("any.fasta")
        .arg("--decompress");
    
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("not yet implemented"));
    
    Ok(())
}