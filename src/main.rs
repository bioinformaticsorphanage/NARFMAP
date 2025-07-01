use std::path::PathBuf;
use std::process;

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use log::{debug, info, warn, error};

mod reference;
mod align;
mod io;
mod hashtable;
mod config;
mod utils;

// Comment out C++ interop temporarily
// extern "C" {
//     fn snappy_max_compressed_length(source_length: size_t) -> size_t;
// }

// #[cxx::bridge]
// pub(crate) mod ffi {
//     unsafe extern "C++" {
//         include!("include/workflow/GenHashTableWorkflow.hpp");
//         include!("include/workflow/Input2SamWorkflow.hpp");
//     }
// }

#[derive(Parser)]
#[command(
    name = "narfmap",
    about = "NARFMAP - A Rust implementation of the Dragen mapper/aligner",
    version,
    author,
    long_about = None
)]
struct Cli {
    /// Set verbosity level (can be used multiple times)
    #[arg(short, long, action = clap::ArgAction::Count)]
    verbose: u8,

    /// Output directory for all generated files
    #[arg(short = 'o', long, value_name = "DIR")]
    output_dir: Option<PathBuf>,
    
    /// Prefix for output files
    #[arg(short = 'p', long, value_name = "PREFIX", default_value = "out")]
    output_prefix: String,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Build a hash table from a reference genome
    BuildHashTable {
        /// Path to the reference FASTA file
        #[arg(short, long, value_name = "FILE", required = true)]
        reference: PathBuf,
        
        /// Hash table k-mer size
        #[arg(short, long, value_name = "INT", default_value = "21")]
        kmer_size: usize,
        
        /// Number of threads to use for hash table generation
        #[arg(short, long, value_name = "INT", default_value = "4")]
        threads: usize,
        
        /// Decompress an existing hash table
        #[arg(long)]
        decompress: bool,
    },

    /// Align reads to a reference genome
    Align {
        /// Path to the reference directory containing the hash table
        #[arg(short, long, value_name = "DIR", required = true)]
        reference_dir: PathBuf,
        
        /// Path to the first FASTQ file with reads
        #[arg(short = '1', long, value_name = "FILE", required = true)]
        fastq1: PathBuf,
        
        /// Path to the second FASTQ file for paired-end reads
        #[arg(short = '2', long, value_name = "FILE")]
        fastq2: Option<PathBuf>,
        
        /// Number of threads to use for alignment
        #[arg(short, long, value_name = "INT", default_value = "4")]
        threads: usize,
        
        /// Read Group ID for SAM output
        #[arg(long, value_name = "ID", default_value = "1")]
        rgid: String,
        
        /// Read Group sample name for SAM output
        #[arg(long, value_name = "SM")]
        rgsm: Option<String>,
    },

    /// Show information about genomic files (FASTA/FASTQ)
    Info {
        /// Path to a FASTA or FASTQ file
        #[arg(short, long, value_name = "FILE", required = true)]
        input: PathBuf,
        
        /// Show detailed information about the first N sequences
        #[arg(short, long, value_name = "INT", default_value = "5")]
        num_sequences: usize,
    },
}

fn setup_logging(verbosity: u8) {
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

fn run_build_hash_table(args: &Commands) -> Result<()> {
    if let Commands::BuildHashTable { reference, kmer_size, threads, decompress } = args {
        info!("Building hash table from reference: {}", reference.display());
        info!("Using k-mer size: {}", kmer_size);
        info!("Using {} threads", threads);
        
        if *decompress {
            info!("Decompressing existing hash table");
            // TODO: Implement hash table decompression
            return Err(anyhow::anyhow!("Hash table decompression not yet implemented"));
        }
        
        // Create hash table configuration
        let mut hash_config = config::HashTableConfig::default();
        hash_config.seed_len = *kmer_size;
        hash_config.num_threads = *threads;
        
        // Create output directory next to reference file
        let output_dir = reference.parent()
            .unwrap_or_else(|| std::path::Path::new("."))
            .join("hash_table");
        
        info!("Output directory: {}", output_dir.display());
        
        // Build the hash table
        let builder = hashtable::HashTableBuilder::new(hash_config)
            .with_context(|| "Failed to create hash table builder")?;
        
        let _hash_table = builder.build_from_fasta(reference, &output_dir)
            .with_context(|| format!("Failed to build hash table from reference: {}", reference.display()))?;
        
        info!("Hash table generation complete");
        info!("Hash table files saved in: {}", output_dir.display());
        Ok(())
    } else {
        unreachable!("Command dispatch error");
    }
}

fn run_align(args: &Commands, output_dir: &Option<PathBuf>, output_prefix: &str) -> Result<()> {
    if let Commands::Align { reference_dir, fastq1, fastq2, threads, rgid, rgsm } = args {
        info!("Aligning reads to reference in directory: {}", reference_dir.display());
        info!("FASTQ file 1: {}", fastq1.display());
        
        // Create alignment configuration
        let mut align_config = config::AlignmentConfig::default();
        align_config.num_threads = *threads;
        
        // Initialize the aligner
        let aligner = align::Aligner::new(align_config, reference_dir)
            .with_context(|| format!("Failed to initialize aligner with reference directory: {}", reference_dir.display()))?;
        
        // Determine output SAM file path
        let output_file = if let Some(dir) = output_dir {
            dir.join(format!("{}.sam", output_prefix))
        } else {
            PathBuf::from(format!("{}.sam", output_prefix))
        };
        
        info!("Output file: {}", output_file.display());
        
        // Create SAM writer
        let mut sam_writer = align::sam::SamWriter::new(&output_file)
            .with_context(|| format!("Failed to create SAM output file: {}", output_file.display()))?;
        
        // Write SAM header (placeholder reference sequences)
        let ref_sequences = vec![("reference".to_string(), 1000000)]; // TODO: Get actual reference info
        sam_writer.write_header(&ref_sequences)?;
        
        if let Some(fastq2) = fastq2 {
            info!("FASTQ file 2: {}", fastq2.display());
            info!("Running paired-end alignment");
            
            // Read paired-end files
            let mut reader1 = io::fastq::FastqReader::from_path(fastq1)
                .with_context(|| format!("Failed to open FASTQ file: {}", fastq1.display()))?;
            let mut reader2 = io::fastq::FastqReader::from_path(fastq2)
                .with_context(|| format!("Failed to open FASTQ file: {}", fastq2.display()))?;
            
            let mut read_count = 0;
            let mut aligned_count = 0;
            
            // Process reads in pairs
            let iter1 = reader1.iter_sequences();
            let iter2 = reader2.iter_sequences();
            
            for (read1_result, read2_result) in iter1.zip(iter2) {
                let read1 = read1_result?;
                let read2 = read2_result?;
                
                read_count += 2;
                
                // Align both reads
                let (alignment1, alignment2) = aligner.align_paired_reads(&read1, &read2)?;
                
                // Write alignments or unmapped reads
                match alignment1 {
                    Some(ref align) => {
                        sam_writer.write_alignment(align, &read1.sequence_string(), &read1.quality_string())?;
                        aligned_count += 1;
                    }
                    None => {
                        sam_writer.write_unmapped(&read1.id, &read1.sequence_string(), &read1.quality_string())?;
                    }
                }
                
                match alignment2 {
                    Some(ref align) => {
                        sam_writer.write_alignment(align, &read2.sequence_string(), &read2.quality_string())?;
                        aligned_count += 1;
                    }
                    None => {
                        sam_writer.write_unmapped(&read2.id, &read2.sequence_string(), &read2.quality_string())?;
                    }
                }
                
                if read_count % 10000 == 0 {
                    info!("Processed {} reads, {} aligned", read_count, aligned_count);
                }
            }
            
            info!("Processed {} total reads, {} aligned ({:.1}%)", 
                  read_count, aligned_count, (aligned_count as f64 / read_count as f64) * 100.0);
        } else {
            info!("Running single-end alignment");
            
            // Read single-end file
            let mut reader = io::fastq::FastqReader::from_path(fastq1)
                .with_context(|| format!("Failed to open FASTQ file: {}", fastq1.display()))?;
            
            let mut read_count = 0;
            let mut aligned_count = 0;
            
            for read_result in reader.iter_sequences() {
                let read = read_result?;
                read_count += 1;
                
                // Align the read
                let alignment = aligner.align_read(&read)?;
                
                // Write alignment or unmapped read
                match alignment {
                    Some(ref align) => {
                        sam_writer.write_alignment(align, &read.sequence_string(), &read.quality_string())?;
                        aligned_count += 1;
                    }
                    None => {
                        sam_writer.write_unmapped(&read.id, &read.sequence_string(), &read.quality_string())?;
                    }
                }
                
                if read_count % 10000 == 0 {
                    info!("Processed {} reads, {} aligned", read_count, aligned_count);
                }
            }
            
            info!("Processed {} total reads, {} aligned ({:.1}%)", 
                  read_count, aligned_count, (aligned_count as f64 / read_count as f64) * 100.0);
        }
        
        sam_writer.flush()?;
        info!("Alignment complete. Output written to: {}", output_file.display());
        Ok(())
    } else {
        unreachable!("Command dispatch error");
    }
}

/// Run the info command to provide information about genomic files
fn run_info(input: &PathBuf, num_sequences: usize) -> Result<()> {
    let path_str = input.to_string_lossy();
    
    if path_str.ends_with(".fasta") || path_str.ends_with(".fa") || path_str.ends_with(".fna") {
        // Handle FASTA file
        println!("Analyzing FASTA file: {}", path_str);
        let sequences = crate::io::fasta::load_reference(input)?;
        
        println!("FASTA file contains {} sequences", sequences.len());
        
        for (i, seq) in sequences.into_iter().take(num_sequences).enumerate() {
            println!("Sequence #{}: {} ({} bp)", i+1, seq.id, seq.len());
        }
    } else if path_str.ends_with(".fastq") || path_str.ends_with(".fq") {
        // Handle FASTQ file
        println!("Analyzing FASTQ file: {}", path_str);
        let mut reader = crate::io::fastq::FastqReader::from_path(input)?;
        
        let mut count = 0;
        let mut total_length = 0;
        let mut displayed = 0;
        
        for seq_result in reader.iter_sequences() {
            let seq = seq_result?;
            count += 1;
            total_length += seq.len();
            
            if displayed < num_sequences {
                println!("Read #{}: {} ({} bp)", displayed+1, seq.id, seq.len());
                displayed += 1;
            }
        }
        
        let avg_length = if count > 0 { total_length as f64 / count as f64 } else { 0.0 };
        println!("FASTQ file contains {} reads with average length {:.1} bp", count, avg_length);
    } else {
        return Err(anyhow::anyhow!("Unsupported file format. Please provide a FASTA or FASTQ file."));
    }
    
    Ok(())
}

fn run() -> Result<()> {
    let cli = Cli::parse();
    
    // Setup logging based on verbosity
    setup_logging(cli.verbose);
    
    // Print version information
    info!("NARFMAP v{}", env!("CARGO_PKG_VERSION"));
    
    // Create output directory if specified
    if let Some(dir) = &cli.output_dir {
        if !dir.exists() {
            std::fs::create_dir_all(dir)
                .with_context(|| format!("Failed to create output directory: {}", dir.display()))?;
            info!("Created output directory: {}", dir.display());
        }
    }
    
    // Dispatch to the appropriate command handler
    match &cli.command {
        Commands::BuildHashTable { .. } => run_build_hash_table(&cli.command),
        Commands::Align { .. } => run_align(&cli.command, &cli.output_dir, &cli.output_prefix),
        Commands::Info { input, num_sequences } => {
            run_info(input, *num_sequences)
        },
    }
}

fn main() {
    if let Err(err) = run() {
        error!("Error: {:?}", err);
        process::exit(1);
    }
}
