use std::path::PathBuf;
use clap::{Parser, Subcommand, ValueEnum};
use colored::*;

#[derive(Parser)]
#[command(
    name = "narfmap",
    about = "A modern, user-friendly DNA sequence aligner",
    version,
    author,
    after_help = "EXAMPLES:
    # Index a reference genome
    narfmap index genome.fasta
    
    # Align single-end reads
    narfmap map genome reads.fastq
    
    # Align paired-end reads
    narfmap map genome reads_1.fastq reads_2.fastq
    
    # Show stats about files
    narfmap stats reads.fastq
    
    # Run with custom settings
    narfmap map genome reads.fastq --threads 16 --output aligned.sam
    
    # Use a config file
    narfmap map genome reads.fastq --config my_settings.toml"
)]
pub struct Cli {
    /// Increase logging verbosity (use multiple times for more detail)
    #[arg(short, long, action = clap::ArgAction::Count, global = true)]
    pub verbose: u8,

    /// Suppress all output except errors
    #[arg(short, long, global = true, conflicts_with = "verbose")]
    pub quiet: bool,

    /// Use a configuration file for settings
    #[arg(short, long, value_name = "FILE", global = true)]
    pub config: Option<PathBuf>,

    /// Disable colored output
    #[arg(long, global = true)]
    pub no_color: bool,

    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Build an index from a reference genome
    #[command(visible_alias = "build")]
    Index {
        /// Reference genome file (FASTA format)
        reference: PathBuf,

        /// Output directory for the index (default: reference_file.index/)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// K-mer size for indexing
        #[arg(short, long, default_value = "21", value_parser = validate_kmer_size)]
        kmer_size: usize,

        /// Number of threads to use
        #[arg(short, long, default_value_t = num_cpus::get())]
        threads: usize,

        /// Force rebuild even if index exists
        #[arg(short, long)]
        force: bool,

        /// Show progress bar during indexing
        #[arg(long, default_value = "true")]
        progress: bool,
    },

    /// Align sequences to a reference genome
    Map {
        /// Reference genome or index directory
        reference: PathBuf,

        /// Input FASTQ file(s) - provide two files for paired-end reads
        #[arg(required = true, num_args = 1..=2)]
        reads: Vec<PathBuf>,

        /// Output file (default: stdout)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Output format
        #[arg(short = 'f', long, default_value = "sam")]
        format: OutputFormat,

        /// Number of threads to use
        #[arg(short, long, default_value_t = num_cpus::get())]
        threads: usize,

        /// Minimum alignment score
        #[arg(long, default_value = "30")]
        min_score: i32,

        /// Minimum mapping quality
        #[arg(long, default_value = "10")]
        min_mapq: u8,

        /// Sample name for SAM header
        #[arg(long)]
        sample: Option<String>,

        /// Read group ID
        #[arg(long, default_value = "1")]
        read_group: String,

        /// Show progress during alignment
        #[arg(long, default_value = "true")]
        progress: bool,

        /// Report alignment statistics
        #[arg(long, default_value = "true")]
        stats: bool,

        /// Save detailed statistics to file
        #[arg(long)]
        stats_file: Option<PathBuf>,
    },

    /// Show statistics about sequence files
    Stats {
        /// Input file(s) to analyze
        #[arg(required = true)]
        files: Vec<PathBuf>,

        /// Show detailed statistics
        #[arg(short, long)]
        detailed: bool,

        /// Number of sequences to preview
        #[arg(short, long, default_value = "5")]
        preview: usize,

        /// Show quality score distribution (FASTQ only)
        #[arg(long)]
        quality: bool,

        /// Output format for stats
        #[arg(short = 'f', long, default_value = "text")]
        format: StatsFormat,
    },

    /// Validate and check file integrity
    Check {
        /// Files to check
        #[arg(required = true)]
        files: Vec<PathBuf>,

        /// Check paired-end read consistency
        #[arg(short, long)]
        paired: bool,

        /// Maximum number of errors to report
        #[arg(long, default_value = "10")]
        max_errors: usize,
    },

    /// Interactive configuration wizard
    Config {
        /// Output configuration file
        #[arg(short, long, default_value = "narfmap.toml")]
        output: PathBuf,

        /// Start with defaults for a specific use case
        #[arg(short, long)]
        preset: Option<Preset>,
    },

    /// Show information about an existing index
    #[command(visible_alias = "index-info")]
    Inspect {
        /// Index directory to inspect
        index: PathBuf,

        /// Show detailed information
        #[arg(short, long)]
        detailed: bool,

        /// Verify index integrity
        #[arg(long)]
        verify: bool,
    },
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum OutputFormat {
    /// SAM format (default)
    Sam,
    /// BAM format (compressed)
    Bam,
    /// CRAM format (reference-compressed)
    Cram,
    /// Tab-delimited format
    Tsv,
    /// PAF format (minimap2-style)
    Paf,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum StatsFormat {
    /// Human-readable text
    Text,
    /// JSON format
    Json,
    /// CSV format
    Csv,
    /// Markdown table
    Markdown,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum Preset {
    /// Short reads (Illumina)
    ShortReads,
    /// Long reads (PacBio/ONT)
    LongReads,
    /// RNA-seq alignment
    RnaSeq,
    /// Bisulfite sequencing
    Bisulfite,
    /// Ancient DNA
    AncientDna,
    /// Metagenomics
    Metagenomics,
}

/// Validate k-mer size is reasonable
fn validate_kmer_size(s: &str) -> Result<usize, String> {
    let size: usize = s.parse()
        .map_err(|_| format!("'{}' is not a valid number", s))?;
    
    if size < 11 {
        Err("K-mer size must be at least 11".to_string())
    } else if size > 31 {
        Err("K-mer size must be at most 31".to_string())
    } else if size % 2 == 0 {
        Err("K-mer size must be odd".to_string())
    } else {
        Ok(size)
    }
}

/// Configuration structure for TOML files
#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct Config {
    pub alignment: AlignmentConfig,
    pub output: OutputConfig,
    pub performance: PerformanceConfig,
}

#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct AlignmentConfig {
    pub min_score: i32,
    pub min_mapq: u8,
    pub max_mismatches: u32,
    pub max_gap_opens: u32,
    pub seed_length: usize,
    pub max_seed_hits: usize,
}

#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct OutputConfig {
    pub format: String,
    pub compression_level: u8,
    pub include_unmapped: bool,
    pub include_secondary: bool,
    pub sort_by_name: bool,
}

#[derive(Debug, serde::Deserialize, serde::Serialize)]
pub struct PerformanceConfig {
    pub threads: usize,
    pub batch_size: usize,
    pub memory_limit: String,
    pub io_threads: usize,
}

impl Default for Config {
    fn default() -> Self {
        Config {
            alignment: AlignmentConfig {
                min_score: 30,
                min_mapq: 10,
                max_mismatches: 5,
                max_gap_opens: 2,
                seed_length: 21,
                max_seed_hits: 500,
            },
            output: OutputConfig {
                format: "sam".to_string(),
                compression_level: 6,
                include_unmapped: true,
                include_secondary: false,
                sort_by_name: false,
            },
            performance: PerformanceConfig {
                threads: num_cpus::get(),
                batch_size: 10000,
                memory_limit: "8G".to_string(),
                io_threads: 2,
            },
        }
    }
}

/// Print a success message with color
pub fn print_success(msg: &str) {
    if atty::is(atty::Stream::Stdout) {
        println!("{} {}", "✓".green().bold(), msg);
    } else {
        println!("✓ {}", msg);
    }
}

/// Print an error message with color
pub fn print_error(msg: &str) {
    if atty::is(atty::Stream::Stderr) {
        eprintln!("{} {}", "✗".red().bold(), msg);
    } else {
        eprintln!("✗ {}", msg);
    }
}

/// Print a warning message with color
pub fn print_warning(msg: &str) {
    if atty::is(atty::Stream::Stderr) {
        eprintln!("{} {}", "⚠".yellow().bold(), msg);
    } else {
        eprintln!("⚠ {}", msg);
    }
}

/// Print an info message with color
pub fn print_info(msg: &str) {
    if atty::is(atty::Stream::Stdout) {
        println!("{} {}", "ℹ".blue().bold(), msg);
    } else {
        println!("ℹ {}", msg);
    }
}

/// Format file size in human-readable form
pub fn format_size(bytes: u64) -> String {
    const UNITS: &[&str] = &["B", "KB", "MB", "GB", "TB"];
    let mut size = bytes as f64;
    let mut unit_idx = 0;
    
    while size >= 1024.0 && unit_idx < UNITS.len() - 1 {
        size /= 1024.0;
        unit_idx += 1;
    }
    
    if unit_idx == 0 {
        format!("{} {}", size as u64, UNITS[unit_idx])
    } else {
        format!("{:.1} {}", size, UNITS[unit_idx])
    }
}

/// Format duration in human-readable form
pub fn format_duration(secs: f64) -> String {
    if secs < 1.0 {
        format!("{:.0}ms", secs * 1000.0)
    } else if secs < 60.0 {
        format!("{:.1}s", secs)
    } else if secs < 3600.0 {
        let mins = secs / 60.0;
        format!("{:.1}m", mins)
    } else {
        let hours = secs / 3600.0;
        format!("{:.1}h", hours)
    }
}