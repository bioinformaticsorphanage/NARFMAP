//! NARFMAP - DRAGEN-OS mapper/aligner rewritten in Rust
//!
//! A high-performance genomic sequence aligner.

use std::ffi::{CStr, CString};
use std::os::raw::c_char;
use std::path::PathBuf;

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};

// Use library modules
use narfmap::alignment::Aligner;
use narfmap::io::{FastqReader, SamWriter};
use narfmap::reference::{HashTable, ReferenceSequence};

// FFI module is only used in main
mod ffi;

// ============================================================================
// CLI Definition
// ============================================================================

#[derive(Parser)]
#[command(name = "narfmap")]
#[command(author = "Illumina/Edico Genome, Rust port")]
#[command(version = env!("CARGO_PKG_VERSION"))]
#[command(about = "DRAGEN-OS mapper/aligner", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,

    /// Verbose output
    #[arg(short, long, global = true)]
    verbose: bool,
}

#[derive(Subcommand)]
enum Commands {
    /// Build hash table from reference FASTA
    #[command(name = "build-hash-table")]
    BuildHashTable {
        /// Reference FASTA file
        #[arg(short = 'r', long = "ht-reference")]
        reference: PathBuf,

        /// Output directory
        #[arg(short = 'o', long = "output-directory")]
        output_dir: PathBuf,

        /// Primary seed length (default: 21)
        #[arg(long = "ht-seed-len", default_value = "21")]
        seed_len: u32,

        /// Maximum seed frequency (default: 16)
        #[arg(long = "ht-max-seed-freq", default_value = "16")]
        max_seed_freq: u32,

        /// Number of threads (default: 8)
        #[arg(long = "ht-num-threads", default_value = "8")]
        num_threads: i32,

        /// Hash table size (e.g., "32GB")
        #[arg(long = "ht-size")]
        size: Option<String>,

        /// Memory limit (e.g., "32GB")
        #[arg(long = "ht-mem-limit")]
        mem_limit: Option<String>,

        /// Decoy sequences FASTA
        #[arg(long = "ht-decoys")]
        decoys: Option<PathBuf>,

        /// BED file for masking regions
        #[arg(long = "ht-mask-bed")]
        mask_bed: Option<PathBuf>,

        /// ALT liftover SAM file
        #[arg(long = "ht-alt-liftover")]
        alt_liftover: Option<PathBuf>,
    },

    /// Align reads to reference
    Align {
        /// Reference directory (containing hash table)
        #[arg(short = 'r', long = "ref-dir")]
        ref_dir: PathBuf,

        /// Input FASTQ file 1
        #[arg(short = '1', long = "fastq1")]
        fastq1: PathBuf,

        /// Input FASTQ file 2 (for paired-end)
        #[arg(short = '2', long = "fastq2")]
        fastq2: Option<PathBuf>,

        /// Output directory
        #[arg(short = 'o', long = "output-directory")]
        output_dir: Option<PathBuf>,

        /// Output file prefix
        #[arg(long = "output-file-prefix", default_value = "output")]
        output_prefix: String,

        /// Read group ID
        #[arg(long = "RGID", default_value = "1")]
        rgid: String,

        /// Read group sample name
        #[arg(long = "RGSM", default_value = "sample")]
        rgsm: String,

        /// Number of threads
        #[arg(long = "num-threads", default_value = "8")]
        num_threads: usize,

        /// Interleaved FASTQ input
        #[arg(long = "interleaved")]
        interleaved: bool,
    },
}

// ============================================================================
// Main Entry Point
// ============================================================================

fn main() -> Result<()> {
    env_logger::init();

    let cli = Cli::parse();

    match cli.command {
        Commands::BuildHashTable {
            reference,
            output_dir,
            seed_len,
            max_seed_freq,
            num_threads,
            size,
            mem_limit,
            decoys,
            mask_bed,
            alt_liftover,
        } => {
            build_hash_table(
                &reference,
                &output_dir,
                seed_len,
                max_seed_freq,
                num_threads,
                size.as_deref(),
                mem_limit.as_deref(),
                decoys.as_deref(),
                mask_bed.as_deref(),
                alt_liftover.as_deref(),
                cli.verbose,
            )?;
        }
        Commands::Align {
            ref_dir,
            fastq1,
            fastq2,
            output_dir,
            output_prefix,
            rgid,
            rgsm,
            num_threads,
            interleaved,
        } => {
            align_reads(
                &ref_dir,
                &fastq1,
                fastq2.as_deref(),
                output_dir.as_deref(),
                &output_prefix,
                &rgid,
                &rgsm,
                num_threads,
                interleaved,
                cli.verbose,
            )?;
        }
    }

    Ok(())
}

// ============================================================================
// Build Hash Table
// ============================================================================

fn build_hash_table(
    reference: &std::path::Path,
    output_dir: &std::path::Path,
    seed_len: u32,
    max_seed_freq: u32,
    num_threads: i32,
    size: Option<&str>,
    mem_limit: Option<&str>,
    decoys: Option<&std::path::Path>,
    mask_bed: Option<&std::path::Path>,
    alt_liftover: Option<&std::path::Path>,
    verbose: bool,
) -> Result<()> {
    use ffi::hash_gen::*;

    println!("NARFMAP Hash Table Builder v{}", env!("CARGO_PKG_VERSION"));
    println!("Building hash table from: {}", reference.display());
    println!("Output directory: {}", output_dir.display());

    // Ensure output directory exists
    std::fs::create_dir_all(output_dir).with_context(|| {
        format!(
            "Failed to create output directory: {}",
            output_dir.display()
        )
    })?;

    // Canonicalize paths
    let reference = reference
        .canonicalize()
        .with_context(|| format!("Reference file not found: {}", reference.display()))?;
    let output_dir = output_dir.canonicalize()?;

    // Convert paths to C strings
    let ref_cstr = path_to_cstring(&reference)?;
    let output_cstr = path_to_cstring(&output_dir)?;

    unsafe {
        // Allocate and initialize config
        let mut config: HashTableConfig = std::mem::zeroed();
        let hdr = Box::into_raw(Box::new(std::mem::zeroed::<HashTableHeader>()));
        config.hdr = hdr;

        // Set default parameters
        set_default_hash_params(&mut config, output_cstr.as_ptr(), HashTableType::Normal);

        // Set reference input
        config.ref_input = strdup(ref_cstr.as_ptr());

        // Set user parameters
        (*config.hdr).pri_seed_bases = seed_len;
        (*config.hdr).max_seed_freq = max_seed_freq;
        config.max_threads = num_threads;

        // Set optional parameters
        if let Some(size_str) = size {
            let cstr = CString::new(size_str)?;
            config.size_str = strdup(cstr.as_ptr());
        }

        if let Some(mem_str) = mem_limit {
            let cstr = CString::new(mem_str)?;
            config.mem_size_str = strdup(cstr.as_ptr());
        }

        if let Some(decoys_path) = decoys {
            let cstr = path_to_cstring(decoys_path)?;
            config.decoy_fname = strdup(cstr.as_ptr());
        }

        if let Some(mask_path) = mask_bed {
            let cstr = path_to_cstring(mask_path)?;
            config.mask_bed = strdup(cstr.as_ptr());
        }

        if let Some(liftover_path) = alt_liftover {
            let cstr = path_to_cstring(liftover_path)?;
            config.alt_liftover = strdup(cstr.as_ptr());
        }

        config.write_hash_file = 1;
        config.write_comp_file = 0;

        if verbose {
            config.show_int_params = 1;
        }

        // Build command line for logging
        let cmd_line = CString::new("narfmap build-hash-table")?;
        config.cmd_line = strdup(cmd_line.as_ptr());

        // Generate hash table
        println!("Generating hash table...");
        let argv: [*mut c_char; 1] = [std::ptr::null_mut()];
        let error_msg = generate_hash_table(&mut config, 0, argv.as_ptr() as *mut *mut c_char);

        // Check for errors
        if !error_msg.is_null() {
            let err_str = CStr::from_ptr(error_msg).to_string_lossy();
            // Free resources
            free_hash_params(&mut config);
            drop(Box::from_raw(hdr));
            bail!("Hash table generation failed: {}", err_str);
        }

        // Cleanup
        free_hash_params(&mut config);
        drop(Box::from_raw(hdr));
    }

    println!("Hash table built successfully!");
    println!("Output files:");
    println!("  - {}/hash_table.bin", output_dir.display());
    println!("  - {}/hash_table.cfg.bin", output_dir.display());
    println!("  - {}/reference.bin", output_dir.display());

    Ok(())
}

fn path_to_cstring(path: &std::path::Path) -> Result<CString> {
    let s = path
        .to_str()
        .ok_or_else(|| anyhow::anyhow!("Path contains invalid UTF-8: {}", path.display()))?;
    CString::new(s).context("Path contains null byte")
}

// ============================================================================
// Align Reads
// ============================================================================

fn align_reads(
    ref_dir: &std::path::Path,
    fastq1: &std::path::Path,
    fastq2: Option<&std::path::Path>,
    output_dir: Option<&std::path::Path>,
    output_prefix: &str,
    rgid: &str,
    rgsm: &str,
    _num_threads: usize,
    _interleaved: bool,
    verbose: bool,
) -> Result<()> {
    println!("NARFMAP Aligner v{}", env!("CARGO_PKG_VERSION"));
    println!("Reference directory: {}", ref_dir.display());
    println!("Input FASTQ 1: {}", fastq1.display());
    if let Some(fq2) = fastq2 {
        println!("Input FASTQ 2: {}", fq2.display());
    }

    // Verify inputs exist
    if !ref_dir.exists() {
        bail!("Reference directory not found: {}", ref_dir.display());
    }
    if !fastq1.exists() {
        bail!("FASTQ file not found: {}", fastq1.display());
    }
    if let Some(fq2) = fastq2 {
        if !fq2.exists() {
            bail!("FASTQ file not found: {}", fq2.display());
        }
    }

    // Load reference and hash table
    println!("Loading reference...");
    let hash_table = HashTable::load(ref_dir).context("Failed to load hash table")?;
    let reference =
        ReferenceSequence::load(ref_dir).context("Failed to load reference sequence")?;

    if verbose {
        println!("  Hash table version: {}", hash_table.config.version);
        println!("  Seed length: {} bp", hash_table.config.seed_len());
        println!("  Reference length: {} bp", reference.total_length());
    }

    // Determine output path
    let out_dir = output_dir.unwrap_or(std::path::Path::new("."));
    let sam_path = out_dir.join(format!("{}.sam", output_prefix));

    // Create output directory if needed
    if let Some(parent) = sam_path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    // Create SAM writer
    println!("Output: {}", sam_path.display());
    let mut sam_writer = SamWriter::create(&sam_path)?;
    sam_writer.write_header(&reference, rgid, rgsm)?;

    // Create aligner
    let aligner = Aligner::new(&hash_table, &reference);

    // Process reads
    println!("Aligning reads...");
    let mut fastq_reader = FastqReader::open(fastq1)?;

    let mut total_reads = 0u64;
    let mut aligned_reads = 0u64;

    for read_result in fastq_reader {
        let read = read_result?;
        total_reads += 1;

        let alignment = aligner.align(&read);

        if alignment.flag & 4 == 0 {
            // Mapped
            aligned_reads += 1;
        }

        sam_writer.write_alignment(&alignment, rgid)?;

        // Progress
        if total_reads % 10000 == 0 {
            eprint!("\r  Processed {} reads...", total_reads);
        }
    }

    eprintln!("\r  Processed {} reads    ", total_reads);

    sam_writer.finish()?;

    // Print summary
    let align_rate = if total_reads > 0 {
        (aligned_reads as f64 / total_reads as f64) * 100.0
    } else {
        0.0
    };

    println!("\nAlignment complete!");
    println!("  Total reads:   {}", total_reads);
    println!("  Aligned reads: {} ({:.1}%)", aligned_reads, align_rate);
    println!("  Output: {}", sam_path.display());

    Ok(())
}
