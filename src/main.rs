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
        } else {
            // TODO: Implement hash table building
        }
        
        info!("Hash table generation complete");
        Ok(())
    } else {
        unreachable!("Command dispatch error");
    }
}

fn run_align(args: &Commands, output_dir: &Option<PathBuf>, output_prefix: &str) -> Result<()> {
    if let Commands::Align { reference_dir, fastq1, fastq2, threads, rgid, rgsm } = args {
        info!("Aligning reads to reference in directory: {}", reference_dir.display());
        info!("FASTQ file 1: {}", fastq1.display());
        
        if let Some(fastq2) = fastq2 {
            info!("FASTQ file 2: {}", fastq2.display());
            info!("Running paired-end alignment");
            // TODO: Implement paired-end alignment
        } else {
            info!("Running single-end alignment");
            // TODO: Implement single-end alignment
        }
        
        info!("Using {} threads", threads);
        info!("Read Group ID: {}", rgid);
        if let Some(sample) = rgsm {
            info!("Read Group Sample: {}", sample);
        }
        
        // TODO: Implement alignment logic
        
        info!("Alignment complete");
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
        info!("Analyzing FASTA file: {}", path_str);
        let sequences = crate::io::fasta::load_reference(input)?;
        
        info!("FASTA file contains {} sequences", sequences.len());
        
        for (i, seq) in sequences.into_iter().take(num_sequences).enumerate() {
            info!("Sequence #{}: {} ({} bp)", i+1, seq.id, seq.len());
        }
    } else if path_str.ends_with(".fastq") || path_str.ends_with(".fq") {
        // Handle FASTQ file
        info!("Analyzing FASTQ file: {}", path_str);
        let mut reader = crate::io::fastq::FastqReader::from_path(input)?;
        
        let mut count = 0;
        let mut total_length = 0;
        let mut displayed = 0;
        
        for seq_result in reader.iter_sequences() {
            let seq = seq_result?;
            count += 1;
            total_length += seq.len();
            
            if displayed < num_sequences {
                info!("Read #{}: {} ({} bp)", displayed+1, seq.id, seq.len());
                displayed += 1;
            }
        }
        
        let avg_length = if count > 0 { total_length as f64 / count as f64 } else { 0.0 };
        info!("FASTQ file contains {} reads with average length {:.1} bp", count, avg_length);
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
