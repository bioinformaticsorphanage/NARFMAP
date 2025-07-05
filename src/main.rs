use std::path::{Path, PathBuf};
use std::process;
use std::time::Instant;

use anyhow::{Context, Result};
use clap::Parser;
use colored::*;
use indicatif::{ProgressBar, ProgressStyle, MultiProgress};
use log::{debug, info, error, warn};

mod cli;
use cli::{Cli, Commands, OutputFormat, StatsFormat, Preset, Config as CliConfig, 
          print_success, print_error, print_warning, print_info, format_size, format_duration};

mod reference;
mod align;
mod io;
mod hashtable;
mod config;
mod utils;

fn main() {
    // Set up color support
    colored::control::set_override(!atty::is(atty::Stream::Stdout));
    
    if let Err(err) = run() {
        print_error(&format!("{:?}", err));
        process::exit(1);
    }
}

fn run() -> Result<()> {
    let cli = Cli::parse();
    
    // Disable colors if requested
    if cli.no_color {
        colored::control::set_override(false);
    }
    
    // Setup logging
    setup_logging(cli.verbose, cli.quiet);
    
    // Load config file if provided
    let _config = if let Some(config_path) = &cli.config {
        let config_str = std::fs::read_to_string(config_path)
            .with_context(|| format!("Failed to read config file: {}", config_path.display()))?;
        toml::from_str::<CliConfig>(&config_str)
            .with_context(|| format!("Failed to parse config file: {}", config_path.display()))?
    } else {
        CliConfig::default()
    };
    
    // Dispatch commands
    match &cli.command {
        Commands::Index { reference, output, kmer_size, threads, force, progress } => {
            run_index(reference, output.as_deref(), *kmer_size, *threads, *force, *progress)
        }
        Commands::Map { reference, reads, output, format, threads, min_score, min_mapq, 
                       sample, read_group, progress, stats, stats_file } => {
            run_map(reference, reads, output.as_deref(), format, *threads, *min_score, 
                   *min_mapq, sample.as_deref(), read_group, *progress, *stats, stats_file.as_deref())
        }
        Commands::Stats { files, detailed, preview, quality, format } => {
            run_stats(files, *detailed, *preview, *quality, *format)
        }
        Commands::Check { files, paired, max_errors } => {
            run_check(files, *paired, *max_errors)
        }
        Commands::Config { output, preset } => {
            run_config_wizard(output, *preset)
        }
        Commands::Inspect { index, detailed, verify } => {
            run_inspect(index, *detailed, *verify)
        }
    }
}

fn setup_logging(verbosity: u8, quiet: bool) {
    if quiet {
        env_logger::Builder::new()
            .filter_level(log::LevelFilter::Error)
            .format_timestamp(Some(env_logger::fmt::TimestampPrecision::Millis))
            .init();
        return;
    }
    
    let log_level = match verbosity {
        0 => log::LevelFilter::Warn,
        1 => log::LevelFilter::Info,
        2 => log::LevelFilter::Debug,
        _ => log::LevelFilter::Trace,
    };
    
    env_logger::Builder::new()
        .filter_level(log_level)
        .format_timestamp(Some(env_logger::fmt::TimestampPrecision::Millis))
        .init();
    
    debug!("Log level set to {:?}", log_level);
}

fn run_index(reference: &PathBuf, output: Option<&Path>, kmer_size: usize, 
             threads: usize, force: bool, show_progress: bool) -> Result<()> {
    let start_time = Instant::now();
    
    // Determine output directory
    let output_dir = if let Some(dir) = output {
        dir.to_path_buf()
    } else {
        // Default: reference_file.index/
        let mut dir = reference.clone();
        let filename = reference.file_stem()
            .ok_or_else(|| anyhow::anyhow!("Invalid reference filename"))?;
        dir.set_file_name(format!("{}.index", filename.to_string_lossy()));
        dir
    };
    
    // Check if index already exists
    if output_dir.exists() && !force {
        print_warning(&format!("Index already exists at: {}", output_dir.display()));
        print_info("Use --force to rebuild the index");
        return Ok(());
    }
    
    print_info(&format!("Building index for: {}", reference.display()));
    print_info(&format!("Output directory: {}", output_dir.display()));
    print_info(&format!("K-mer size: {}", kmer_size));
    print_info(&format!("Threads: {}", threads));
    
    // Create progress bar
    let progress_bar = if show_progress {
        let pb = ProgressBar::new_spinner();
        pb.set_style(ProgressStyle::default_spinner()
            .template("{spinner:.green} {msg}")
            .unwrap());
        pb.set_message("Loading reference genome...");
        Some(pb)
    } else {
        None
    };
    
    // Create hash table configuration
    let mut hash_config = config::HashTableConfig::default();
    hash_config.seed_len = kmer_size;
    hash_config.num_threads = threads;
    
    // Build the hash table
    let builder = hashtable::HashTableBuilder::new(hash_config)
        .with_context(|| "Failed to create hash table builder")?;
    
    if let Some(pb) = &progress_bar {
        pb.set_message("Building hash table...");
    }
    
    let _hash_table = builder.build_from_fasta(reference, &output_dir)
        .with_context(|| format!("Failed to build hash table from reference: {}", reference.display()))?;
    
    if let Some(pb) = progress_bar {
        pb.finish_and_clear();
    }
    
    let elapsed = start_time.elapsed();
    print_success(&format!("Index built successfully in {}", format_duration(elapsed.as_secs_f64())));
    print_info(&format!("Index location: {}", output_dir.display()));
    
    // Show index size
    if let Ok(size) = get_directory_size(&output_dir) {
        print_info(&format!("Index size: {}", format_size(size)));
    }
    
    Ok(())
}

fn run_map(reference: &PathBuf, reads: &[PathBuf], output: Option<&Path>, 
           _format: &OutputFormat, threads: usize, min_score: i32, min_mapq: u8,
           sample: Option<&str>, read_group: &str, show_progress: bool, 
           show_stats: bool, stats_file: Option<&Path>) -> Result<()> {
    let start_time = Instant::now();
    
    // Determine if this is paired-end
    let is_paired = reads.len() == 2;
    
    print_info(&format!("Aligning {} reads to: {}", 
                       if is_paired { "paired-end" } else { "single-end" },
                       reference.display()));
    
    // Find the index directory
    let index_dir = if reference.is_dir() {
        reference.clone()
    } else {
        // Check for default index location
        let mut dir = reference.clone();
        let filename = reference.file_stem()
            .ok_or_else(|| anyhow::anyhow!("Invalid reference filename"))?;
        dir.set_file_name(format!("{}.index", filename.to_string_lossy()));
        if dir.exists() {
            print_info(&format!("Using index at: {}", dir.display()));
            dir
        } else {
            return Err(anyhow::anyhow!(
                "No index found. Please run 'narfmap index' first or specify the index directory."
            ));
        }
    };
    
    // Create alignment configuration
    let mut align_config = config::AlignmentConfig::default();
    align_config.num_threads = threads;
    align_config.min_score = min_score;
    // Store min_mapq for later use
    let _min_mapq = min_mapq;
    
    // Initialize the aligner
    let mut aligner = align::Aligner::new(align_config, &index_dir)
        .with_context(|| format!("Failed to initialize aligner with index: {}", index_dir.display()))?;
    
    // Start statistics tracking
    aligner.start_stats();
    
    // Create output writer
    let mut sam_writer = if let Some(output_path) = output {
        print_info(&format!("Output file: {}", output_path.display()));
        align::sam::SamWriter::new(output_path)?
    } else {
        print_info("Output: stdout");
        align::sam::SamWriter::new_stdout()?
    };
    
    // Create progress tracking
    let multi_progress = MultiProgress::new();
    let progress_bar = if show_progress {
        let pb = multi_progress.add(ProgressBar::new(0));
        pb.set_style(ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} reads ({per_sec})")
            .unwrap()
            .progress_chars("#>-"));
        Some(pb)
    } else {
        None
    };
    
    // Get reference sequences for header (placeholder for now)
    let ref_sequences = vec![("reference".to_string(), 1000000)]; // TODO: Get actual reference info
    sam_writer.write_header(&ref_sequences)?;
    
    let mut total_reads = 0u64;
    let mut aligned_reads = 0u64;
    
    if is_paired {
        // Paired-end alignment
        print_info(&format!("Read 1: {}", reads[0].display()));
        print_info(&format!("Read 2: {}", reads[1].display()));
        
        let mut reader1 = io::fastq::FastqReader::from_path(&reads[0])?;
        let mut reader2 = io::fastq::FastqReader::from_path(&reads[1])?;
        
        // Count reads for progress bar
        if let Some(pb) = &progress_bar {
            if let Ok(count) = count_fastq_reads(&reads[0]) {
                pb.set_length(count * 2);
            }
        }
        
        let iter1 = reader1.iter_sequences();
        let iter2 = reader2.iter_sequences();
        
        for (read1_result, read2_result) in iter1.zip(iter2) {
            let read1 = read1_result?;
            let read2 = read2_result?;
            
            total_reads += 2;
            
            // Align both reads
            let (alignment1, alignment2) = aligner.align_paired_reads(&read1, &read2)?;
            
            // Write alignments
            match alignment1 {
                Some(ref align) => {
                    sam_writer.write_alignment(align, &read1.sequence_string(), &read1.quality_string())?;
                    aligned_reads += 1;
                }
                None => {
                    sam_writer.write_unmapped(&read1.id, &read1.sequence_string(), &read1.quality_string())?;
                }
            }
            
            match alignment2 {
                Some(ref align) => {
                    sam_writer.write_alignment(align, &read2.sequence_string(), &read2.quality_string())?;
                    aligned_reads += 1;
                }
                None => {
                    sam_writer.write_unmapped(&read2.id, &read2.sequence_string(), &read2.quality_string())?;
                }
            }
            
            if let Some(pb) = &progress_bar {
                pb.set_position(total_reads);
            }
        }
    } else {
        // Single-end alignment
        print_info(&format!("Reads: {}", reads[0].display()));
        
        let mut reader = io::fastq::FastqReader::from_path(&reads[0])?;
        
        // Count reads for progress bar
        if let Some(pb) = &progress_bar {
            if let Ok(count) = count_fastq_reads(&reads[0]) {
                pb.set_length(count);
            }
        }
        
        for read_result in reader.iter_sequences() {
            let read = read_result?;
            total_reads += 1;
            
            // Align the read
            let alignment = aligner.align_read(&read)?;
            
            // Write alignment
            match alignment {
                Some(ref align) => {
                    sam_writer.write_alignment(align, &read.sequence_string(), &read.quality_string())?;
                    aligned_reads += 1;
                }
                None => {
                    sam_writer.write_unmapped(&read.id, &read.sequence_string(), &read.quality_string())?;
                }
            }
            
            if let Some(pb) = &progress_bar {
                pb.set_position(total_reads);
            }
        }
    }
    
    if let Some(pb) = progress_bar {
        pb.finish_and_clear();
    }
    
    sam_writer.flush()?;
    
    // Finalize statistics
    aligner.finalize_stats();
    let stats = aligner.get_stats();
    
    let elapsed = start_time.elapsed();
    let alignment_rate = (aligned_reads as f64 / total_reads as f64) * 100.0;
    
    print_success(&format!("Alignment completed in {}", format_duration(elapsed.as_secs_f64())));
    print_info(&format!("Total reads: {}", total_reads));
    print_info(&format!("Aligned reads: {} ({:.1}%)", aligned_reads, alignment_rate));
    print_info(&format!("Speed: {:.0} reads/sec", total_reads as f64 / elapsed.as_secs_f64()));
    
    if show_stats {
        println!("\n{}", stats.generate_report());
    }
    
    // Save detailed statistics if requested
    if let Some(stats_path) = stats_file {
        if let Ok(json_stats) = stats.to_json() {
            std::fs::write(stats_path, json_stats)
                .with_context(|| format!("Failed to write statistics to: {}", stats_path.display()))?;
            print_info(&format!("Detailed statistics saved to: {}", stats_path.display()));
        }
    }
    
    Ok(())
}

fn run_stats(files: &[PathBuf], detailed: bool, preview: usize, 
             show_quality: bool, format: StatsFormat) -> Result<()> {
    for file in files {
        print_info(&format!("Analyzing: {}", file.display()));
        
        let extension = file.extension()
            .and_then(|e| e.to_str())
            .unwrap_or("");
        
        match extension {
            "fasta" | "fa" | "fna" => analyze_fasta(file, detailed, preview, format)?,
            "fastq" | "fq" => analyze_fastq(file, detailed, preview, show_quality, format)?,
            _ => print_warning(&format!("Unknown file type: {}", file.display())),
        }
        
        if files.len() > 1 {
            println!(); // Separator between files
        }
    }
    
    Ok(())
}

fn analyze_fasta(file: &PathBuf, detailed: bool, preview: usize, format: StatsFormat) -> Result<()> {
    let sequences = io::fasta::load_reference(file)?;
    
    match format {
        StatsFormat::Text => {
            println!("File type: FASTA");
            println!("Sequences: {}", sequences.len());
            
            let total_length: usize = sequences.iter().map(|s| s.len()).sum();
            println!("Total length: {} bp", total_length);
            
            if !sequences.is_empty() {
                let avg_length = total_length / sequences.len();
                let min_length = sequences.iter().map(|s| s.len()).min().unwrap();
                let max_length = sequences.iter().map(|s| s.len()).max().unwrap();
                
                println!("Average length: {} bp", avg_length);
                println!("Min length: {} bp", min_length);
                println!("Max length: {} bp", max_length);
            }
            
            if detailed && !sequences.is_empty() {
                println!("\nSequence details:");
                for (i, seq) in sequences.iter().take(preview).enumerate() {
                    println!("  {}: {} ({} bp)", i + 1, seq.id, seq.len());
                }
                if sequences.len() > preview {
                    println!("  ... and {} more", sequences.len() - preview);
                }
            }
        }
        StatsFormat::Json => {
            let stats = serde_json::json!({
                "file": file.display().to_string(),
                "type": "FASTA",
                "sequences": sequences.len(),
                "total_length": sequences.iter().map(|s| s.len()).sum::<usize>(),
                "sequences_preview": sequences.iter().take(preview).map(|s| {
                    serde_json::json!({
                        "id": s.id,
                        "length": s.len()
                    })
                }).collect::<Vec<_>>()
            });
            println!("{}", serde_json::to_string_pretty(&stats)?);
        }
        _ => print_warning("Format not yet implemented"),
    }
    
    Ok(())
}

fn analyze_fastq(file: &PathBuf, detailed: bool, preview: usize, 
                 show_quality: bool, format: StatsFormat) -> Result<()> {
    let mut reader = io::fastq::FastqReader::from_path(file)?;
    
    let mut count = 0;
    let mut total_length = 0;
    let mut min_length = usize::MAX;
    let mut max_length = 0;
    let mut quality_scores: Vec<u8> = Vec::new();
    let mut preview_seqs = Vec::new();
    
    for seq_result in reader.iter_sequences() {
        let seq = seq_result?;
        count += 1;
        let len = seq.len();
        total_length += len;
        min_length = min_length.min(len);
        max_length = max_length.max(len);
        
        if show_quality {
            // TODO: Implement quality score extraction
            // quality_scores.extend(seq.quality_scores());
        }
        
        if preview_seqs.len() < preview {
            preview_seqs.push((seq.id.clone(), len));
        }
    }
    
    match format {
        StatsFormat::Text => {
            println!("File type: FASTQ");
            println!("Reads: {}", count);
            println!("Total length: {} bp", total_length);
            
            if count > 0 {
                let avg_length = total_length / count;
                println!("Average length: {} bp", avg_length);
                println!("Min length: {} bp", min_length);
                println!("Max length: {} bp", max_length);
            }
            
            if show_quality && !quality_scores.is_empty() {
                // TODO: Fix quality score averaging
                // let avg_qual = quality_scores.iter().map(|&q| q as f64).sum::<f64>() 
                //     / quality_scores.len() as f64;
                // println!("Average quality: {:.1}", avg_qual);
            }
            
            if detailed && !preview_seqs.is_empty() {
                println!("\nRead preview:");
                for (i, (id, len)) in preview_seqs.iter().enumerate() {
                    println!("  {}: {} ({} bp)", i + 1, id, len);
                }
                if count > preview {
                    println!("  ... and {} more", count - preview);
                }
            }
        }
        StatsFormat::Json => {
            let mut stats = serde_json::json!({
                "file": file.display().to_string(),
                "type": "FASTQ",
                "reads": count,
                "total_length": total_length,
                "min_length": min_length,
                "max_length": max_length,
                "average_length": if count > 0 { total_length / count } else { 0 }
            });
            
            if show_quality && !quality_scores.is_empty() {
                // TODO: Fix quality score averaging
                // let avg_qual = quality_scores.iter().map(|&q| q as f64).sum::<f64>() 
                //     / quality_scores.len() as f64;
                // stats["average_quality"] = serde_json::json!(avg_qual);
            }
            
            println!("{}", serde_json::to_string_pretty(&stats)?);
        }
        _ => print_warning("Format not yet implemented"),
    }
    
    Ok(())
}

fn run_check(files: &[PathBuf], check_paired: bool, max_errors: usize) -> Result<()> {
    print_info("Checking file integrity...");
    
    let mut total_errors = 0;
    
    for file in files {
        let mut file_errors = 0;
        print_info(&format!("Checking: {}", file.display()));
        
        let extension = file.extension()
            .and_then(|e| e.to_str())
            .unwrap_or("");
        
        match extension {
            "fastq" | "fq" => {
                let mut reader = io::fastq::FastqReader::from_path(file)?;
                let mut line_num = 0;
                
                for (i, result) in reader.iter_sequences().enumerate() {
                    line_num = i * 4 + 1; // FASTQ has 4 lines per record
                    
                    if let Err(e) = result {
                        file_errors += 1;
                        total_errors += 1;
                        print_error(&format!("Error at line ~{}: {}", line_num, e));
                        
                        if total_errors >= max_errors {
                            print_warning(&format!("Reached maximum error count ({})", max_errors));
                            return Ok(());
                        }
                    }
                }
                
                if file_errors == 0 {
                    print_success(&format!("{}: OK", file.display()));
                } else {
                    print_error(&format!("{}: {} errors found", file.display(), file_errors));
                }
            }
            _ => {
                print_warning(&format!("Check not implemented for file type: {}", extension));
            }
        }
    }
    
    if check_paired && files.len() == 2 {
        print_info("Checking paired-end consistency...");
        // TODO: Implement paired-end consistency check
        print_warning("Paired-end check not yet implemented");
    }
    
    if total_errors == 0 {
        print_success("All files passed validation");
    } else {
        print_error(&format!("Total errors found: {}", total_errors));
    }
    
    Ok(())
}

fn run_config_wizard(output: &PathBuf, preset: Option<Preset>) -> Result<()> {
    print_info("Configuration wizard");
    
    let config = if let Some(preset) = preset {
        print_info(&format!("Using preset: {:?}", preset));
        // TODO: Implement preset configurations
        CliConfig::default()
    } else {
        print_info("Using default configuration");
        CliConfig::default()
    };
    
    let toml_string = toml::to_string_pretty(&config)?;
    std::fs::write(output, toml_string)?;
    
    print_success(&format!("Configuration saved to: {}", output.display()));
    print_info("You can now use this configuration with: --config");
    
    Ok(())
}

fn run_inspect(index: &PathBuf, detailed: bool, verify: bool) -> Result<()> {
    print_info(&format!("Inspecting index: {}", index.display()));
    
    if !index.exists() {
        return Err(anyhow::anyhow!("Index directory not found: {}", index.display()));
    }
    
    if !index.is_dir() {
        return Err(anyhow::anyhow!("Index path is not a directory: {}", index.display()));
    }
    
    // Check for required files
    let config_file = index.join("hash_table.cfg");
    let _bucket_file = index.join("hash_table_buckets_0.dat");
    let _metadata_file = index.join("reference_metadata.json");
    
    if !config_file.exists() {
        return Err(anyhow::anyhow!("Index configuration not found"));
    }
    
    // Load and display configuration
    let config_data = std::fs::read_to_string(&config_file)?;
    println!("Index configuration:");
    println!("{}", config_data);
    
    // Check index size
    if let Ok(size) = get_directory_size(index) {
        print_info(&format!("Index size: {}", format_size(size)));
    }
    
    // List index files
    if detailed {
        println!("\nIndex files:");
        for entry in std::fs::read_dir(index)? {
            let entry = entry?;
            let metadata = entry.metadata()?;
            println!("  {} ({})", 
                     entry.file_name().to_string_lossy(),
                     format_size(metadata.len()));
        }
    }
    
    if verify {
        print_info("Verifying index integrity...");
        // TODO: Implement index verification
        print_warning("Index verification not yet implemented");
    }
    
    print_success("Index appears valid");
    Ok(())
}

// Helper functions

fn count_fastq_reads(file: &PathBuf) -> Result<u64> {
    let mut reader = io::fastq::FastqReader::from_path(file)?;
    let count = reader.iter_sequences().count() as u64;
    Ok(count)
}

fn get_directory_size(path: &PathBuf) -> Result<u64> {
    let mut total_size = 0;
    
    for entry in std::fs::read_dir(path)? {
        let entry = entry?;
        let metadata = entry.metadata()?;
        if metadata.is_file() {
            total_size += metadata.len();
        }
    }
    
    Ok(total_size)
}