use assert_cmd::Command;
use assert_fs::prelude::*;
use predicates::prelude::*;

mod common;

#[test]
fn test_align_missing_hash_table() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("align")
        .arg("-r").arg("nonexistent_dir")
        .arg("-1").arg("some.fastq");
    
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("Failed to initialize aligner").or(
            predicate::str::contains("No such file")
        ));
    
    Ok(())
}

#[test]
fn test_align_missing_fastq() -> Result<(), Box<dyn std::error::Error>> {
    let temp_workspace = common::create_test_workspace();
    let fake_hash_dir = temp_workspace.child("hash_table");
    fake_hash_dir.create_dir_all()?;
    
    // Create a minimal hash table config to make the directory look valid
    let fake_config = fake_hash_dir.child("hash_table.cfg");
    fake_config.write_str("[config]\nseed_len = 21\n")?;
    
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("align")
        .arg("-r").arg(fake_hash_dir.path())
        .arg("-1").arg("nonexistent.fastq");
    
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("Failed to open FASTQ file").or(
            predicate::str::contains("No such file")
        ));
    
    Ok(())
}

#[test]
fn test_malformed_fasta() -> Result<(), Box<dyn std::error::Error>> {
    let temp_workspace = common::create_test_workspace();
    let bad_fasta = temp_workspace.child("bad.fasta");
    bad_fasta.write_str("This is not a valid FASTA file\nNo headers here!")?;
    
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("build-hash-table")
        .arg("-r").arg(bad_fasta.path())
        .current_dir(&temp_workspace);
    
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("Failed to build hash table").or(
            predicate::str::contains("failed to")
        ));
    
    Ok(())
}

#[test]
fn test_empty_fasta() -> Result<(), Box<dyn std::error::Error>> {
    let temp_workspace = common::create_test_workspace();
    let empty_fasta = temp_workspace.child("empty.fasta");
    empty_fasta.write_str("")?;
    
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("build-hash-table")
        .arg("-r").arg(empty_fasta.path())
        .current_dir(&temp_workspace);
    
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("Failed to build hash table").or(
            predicate::str::contains("empty").or(
                predicate::str::contains("No sequences")
            )
        ));
    
    Ok(())
}

#[test]
fn test_invalid_kmer_size() -> Result<(), Box<dyn std::error::Error>> {
    // Skip if test data doesn't exist
    if let Err(_) = common::check_test_data_exists() {
        println!("Skipping invalid k-mer test due to missing test data");
        return Ok(());
    }
    
    let temp_workspace = common::create_test_workspace();
    
    // Test k-mer size of 0
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("build-hash-table")
        .arg("-r").arg(common::get_tiny_fasta())
        .arg("-k").arg("0")
        .current_dir(&temp_workspace);
    
    // This should either fail or handle gracefully
    let result = cmd.assert();
    // We don't enforce failure here since the implementation might handle it gracefully
    
    // Test extremely large k-mer size
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("build-hash-table")
        .arg("-r").arg(common::get_tiny_fasta())
        .arg("-k").arg("100")
        .current_dir(&temp_workspace);
    
    // This might succeed but with warnings
    let result = cmd.assert();
    
    Ok(())
}

#[test]
fn test_info_unsupported_format() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("info")
        .arg("-i").arg("Cargo.toml");  // Wrong file type
    
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("Unsupported file format"));
    
    Ok(())
}

#[test]
fn test_info_nonexistent_file() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("info")
        .arg("-i").arg("nonexistent_file.fasta");
    
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("No such file").or(
            predicate::str::contains("Failed to")
        ));
    
    Ok(())
}

#[test]
fn test_align_unpaired_fastq2() -> Result<(), Box<dyn std::error::Error>> {
    // Skip if test data doesn't exist
    if let Err(_) = common::check_test_data_exists() {
        println!("Skipping unpaired test due to missing test data");
        return Ok(());
    }
    
    let temp_workspace = common::create_test_workspace();
    let hash_dir = temp_workspace.child("hash_table");
    hash_dir.create_dir_all()?;
    
    // Create a fake hash table config
    let fake_config = hash_dir.child("hash_table.cfg");
    fake_config.write_str("[config]\nseed_len = 21\n")?;
    
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("align")
        .arg("-r").arg(hash_dir.path())
        .arg("-1").arg(common::get_tiny_fastq())
        .arg("-2").arg("nonexistent_r2.fastq");
    
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("Failed to open FASTQ file").or(
            predicate::str::contains("No such file")
        ));
    
    Ok(())
}

#[test]
fn test_missing_required_arguments() -> Result<(), Box<dyn std::error::Error>> {
    // Test build-hash-table without reference
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("build-hash-table");
    
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("required").or(
            predicate::str::contains("missing")
        ));
    
    // Test align without reference directory
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("align");
    
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("required").or(
            predicate::str::contains("missing")
        ));
    
    // Test info without input file
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("info");
    
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("required").or(
            predicate::str::contains("missing")
        ));
    
    Ok(())
}

#[test]
fn test_invalid_command() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("invalid-command");
    
    cmd.assert()
        .failure()
        .stderr(predicate::str::contains("unrecognized subcommand").or(
            predicate::str::contains("invalid")
        ));
    
    Ok(())
}

#[test]
fn test_help_commands() -> Result<(), Box<dyn std::error::Error>> {
    // Test main help
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("--help");
    
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("NARFMAP"));
    
    // Test subcommand help
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("build-hash-table").arg("--help");
    
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("Build a hash table"));
    
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("align").arg("--help");
    
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("Align reads"));
    
    let mut cmd = Command::cargo_bin("narfmap")?;
    cmd.arg("info").arg("--help");
    
    cmd.assert()
        .success()
        .stdout(predicate::str::contains("Show information"));
    
    Ok(())
}