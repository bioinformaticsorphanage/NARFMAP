//! NARFMAP - DRAGEN-OS mapper/aligner rewritten in Rust
//!
//! A high-performance genomic sequence aligner.

use std::ffi::{CStr, CString};
use std::os::raw::c_char;
use std::path::PathBuf;

use anyhow::{bail, Context, Result};
use clap::{Parser, Subcommand};
use rayon::prelude::*;

// Use library modules
use narfmap::alignment::Aligner;
use narfmap::io::{FastqReader, SamWriter};
use narfmap::reference::{HashTable, ReferenceSequence};
use narfmap::{Alignment, Read};

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

        /// Minimum insert size for proper pairing
        #[arg(long = "min-insert-size", default_value = "0")]
        min_insert_size: u64,

        /// Maximum insert size for proper pairing
        #[arg(long = "max-insert-size", default_value = "1000")]
        max_insert_size: u64,

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
            min_insert_size,
            max_insert_size,
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
                min_insert_size,
                max_insert_size,
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

#[allow(clippy::too_many_arguments)]
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
    use ffi::hash_gen::{
        free_hash_params, generate_hash_table, set_default_hash_params, strdup, HashTableConfig,
        HashTableHeader, HashTableType,
    };

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

    // SAFETY: FFI calls to C hash table generation library
    #[allow(unsafe_code)]
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
            bail!("Hash table generation failed: {err_str}");
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

#[allow(clippy::too_many_arguments)]
fn align_reads(
    ref_dir: &std::path::Path,
    fastq1: &std::path::Path,
    fastq2: Option<&std::path::Path>,
    output_dir: Option<&std::path::Path>,
    output_prefix: &str,
    rgid: &str,
    rgsm: &str,
    num_threads: usize,
    min_insert_size: u64,
    max_insert_size: u64,
    interleaved: bool,
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
    if interleaved && fastq2.is_some() {
        bail!("--interleaved cannot be used with --fastq2");
    }
    if let Some(fq2) = fastq2 {
        if !fq2.exists() {
            bail!("FASTQ file not found: {}", fq2.display());
        }
    }

    if min_insert_size > max_insert_size {
        bail!("--min-insert-size cannot be greater than --max-insert-size");
    }

    let pairing = PairingHeuristic {
        min_insert_size,
        max_insert_size,
    };

    let thread_count = num_threads.max(1);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build()?;
    let chunk_size = 10_000usize;

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
    let sam_path = out_dir.join(format!("{output_prefix}.sam"));

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
    let mut total_reads = 0u64;
    let mut aligned_reads = 0u64;

    if interleaved {
        let mut fastq_reader = FastqReader::open(fastq1)?;
        loop {
            let pairs = read_chunk_pairs_interleaved(&mut fastq_reader, chunk_size)?;
            if pairs.is_empty() {
                break;
            }

            let results = pool.install(|| {
                pairs
                    .par_iter()
                    .map(|pair| align_pair(&aligner, pair, &pairing))
                    .collect::<Vec<_>>()
            });

            total_reads += (pairs.len() as u64) * 2;
            for (aln1, aln2) in results {
                if aln1.flag & 4 == 0 {
                    aligned_reads += 1;
                }
                if aln2.flag & 4 == 0 {
                    aligned_reads += 1;
                }

                sam_writer.write_alignment(&aln1, rgid)?;
                sam_writer.write_alignment(&aln2, rgid)?;
            }

            if total_reads % 10000 == 0 {
                eprint!("\r  Processed {total_reads} reads...");
            }
        }
    } else if let Some(fq2) = fastq2 {
        let mut fastq_reader1 = FastqReader::open(fastq1)?;
        let mut fastq_reader2 = FastqReader::open(fq2)?;
        loop {
            let pairs = read_chunk_pairs(&mut fastq_reader1, &mut fastq_reader2, chunk_size)?;
            if pairs.is_empty() {
                break;
            }

            let results = pool.install(|| {
                pairs
                    .par_iter()
                    .map(|pair| align_pair(&aligner, pair, &pairing))
                    .collect::<Vec<_>>()
            });

            total_reads += (pairs.len() as u64) * 2;
            for (aln1, aln2) in results {
                if aln1.flag & 4 == 0 {
                    aligned_reads += 1;
                }
                if aln2.flag & 4 == 0 {
                    aligned_reads += 1;
                }

                sam_writer.write_alignment(&aln1, rgid)?;
                sam_writer.write_alignment(&aln2, rgid)?;
            }

            if total_reads % 10000 == 0 {
                eprint!("\r  Processed {total_reads} reads...");
            }
        }
    } else {
        let mut fastq_reader = FastqReader::open(fastq1)?;
        loop {
            let reads = read_chunk(&mut fastq_reader, chunk_size)?;
            if reads.is_empty() {
                break;
            }

            let alignments = pool.install(|| {
                reads
                    .par_iter()
                    .map(|read| aligner.align(read))
                    .collect::<Vec<_>>()
            });

            total_reads += reads.len() as u64;
            for alignment in alignments {
                if alignment.flag & 4 == 0 {
                    aligned_reads += 1;
                }
                sam_writer.write_alignment(&alignment, rgid)?;
            }

            if total_reads % 10000 == 0 {
                eprint!("\r  Processed {total_reads} reads...");
            }
        }
    }

    eprintln!("\r  Processed {total_reads} reads    ");

    sam_writer.finish()?;

    // Print summary
    let align_rate = if total_reads > 0 {
        (aligned_reads as f64 / total_reads as f64) * 100.0
    } else {
        0.0
    };

    println!("\nAlignment complete!");
    println!("  Total reads:   {total_reads}");
    println!("  Aligned reads: {aligned_reads} ({align_rate:.1}%)");
    println!("  Output: {}", sam_path.display());

    Ok(())
}

fn read_chunk(reader: &mut FastqReader, chunk_size: usize) -> Result<Vec<Read>> {
    let mut reads = Vec::with_capacity(chunk_size);
    for _ in 0..chunk_size {
        match reader.read_record()? {
            Some(read) => reads.push(read),
            None => break,
        }
    }
    Ok(reads)
}

fn read_chunk_pairs_interleaved(
    reader: &mut FastqReader,
    chunk_size: usize,
) -> Result<Vec<ReadPair>> {
    let mut pairs = Vec::with_capacity(chunk_size);
    for _ in 0..chunk_size {
        let read1 = match reader.read_record()? {
            Some(read) => read,
            None => break,
        };
        let read2 = reader
            .read_record()?
            .ok_or_else(|| anyhow::anyhow!("Interleaved FASTQ missing mate for {}", read1.name))?;
        let pair_name = pair_name(&read1, &read2)?;
        pairs.push(ReadPair {
            read1,
            read2,
            pair_name,
        });
    }
    Ok(pairs)
}

fn read_chunk_pairs(
    reader1: &mut FastqReader,
    reader2: &mut FastqReader,
    chunk_size: usize,
) -> Result<Vec<ReadPair>> {
    let mut pairs = Vec::with_capacity(chunk_size);
    for _ in 0..chunk_size {
        let read1 = reader1.read_record()?;
        let read2 = reader2.read_record()?;
        match (read1, read2) {
            (None, None) => break,
            (Some(_), None) => {
                bail!("FASTQ 2 ended before FASTQ 1")
            }
            (None, Some(_)) => {
                bail!("FASTQ 1 ended before FASTQ 2")
            }
            (Some(read1), Some(read2)) => {
                let pair_name = pair_name(&read1, &read2)?;
                pairs.push(ReadPair {
                    read1,
                    read2,
                    pair_name,
                });
            }
        }
    }
    Ok(pairs)
}

fn pair_name(read1: &Read, read2: &Read) -> Result<String> {
    let name1 = normalize_pair_name(&read1.name);
    let name2 = normalize_pair_name(&read2.name);
    if name1 != name2 {
        bail!("FASTQ pair name mismatch: {} vs {}", read1.name, read2.name);
    }
    Ok(name1)
}

fn normalize_pair_name(name: &str) -> String {
    if let Some(stripped) = name.strip_suffix("/1") {
        stripped.to_string()
    } else if let Some(stripped) = name.strip_suffix("/2") {
        stripped.to_string()
    } else {
        name.to_string()
    }
}

#[derive(Clone, Copy, Debug)]
struct PairingHeuristic {
    min_insert_size: u64,
    max_insert_size: u64,
}

impl PairingHeuristic {
    fn is_proper_pair(&self, aln1: &Alignment, aln2: &Alignment, template_length: i64) -> bool {
        if template_length == 0 {
            return false;
        }
        if aln1.flag & 4 != 0 || aln2.flag & 4 != 0 {
            return false;
        }
        if aln1.ref_name != aln2.ref_name {
            return false;
        }

        let tlen_abs = template_length.unsigned_abs();
        if tlen_abs < self.min_insert_size || tlen_abs > self.max_insert_size {
            return false;
        }

        let read1_forward = aln1.flag & 0x10 == 0;
        let read2_reverse = aln2.flag & 0x10 != 0;
        if !read1_forward || !read2_reverse {
            return false;
        }

        aln1.position <= aln2.position
    }
}

struct ReadPair {
    read1: Read,
    read2: Read,
    pair_name: String,
}

fn align_pair(
    aligner: &Aligner,
    pair: &ReadPair,
    pairing: &PairingHeuristic,
) -> (Alignment, Alignment) {
    let mut aln1 = aligner.align(&pair.read1);
    let mut aln2 = aligner.align(&pair.read2);
    apply_pair_metadata(
        &mut aln1,
        &mut aln2,
        &pair.read1,
        &pair.read2,
        &pair.pair_name,
        pairing,
    );
    (aln1, aln2)
}

fn apply_pair_metadata(
    aln1: &mut Alignment,
    aln2: &mut Alignment,
    read1: &Read,
    read2: &Read,
    pair_name: &str,
    pairing: &PairingHeuristic,
) {
    let unmapped1 = aln1.flag & 4 != 0;
    let unmapped2 = aln2.flag & 4 != 0;

    aln1.read_name = pair_name.to_string();
    aln2.read_name = pair_name.to_string();

    aln1.flag |= 0x1 | 0x40;
    aln2.flag |= 0x1 | 0x80;

    if unmapped2 {
        aln1.flag |= 0x8;
    }
    if unmapped1 {
        aln2.flag |= 0x8;
    }

    if !unmapped2 && (aln2.flag & 0x10 != 0) {
        aln1.flag |= 0x20;
    }
    if !unmapped1 && (aln1.flag & 0x10 != 0) {
        aln2.flag |= 0x20;
    }

    if unmapped2 {
        aln1.mate_ref_name = "*".to_string();
        aln1.mate_position = 0;
    } else {
        aln1.mate_ref_name = if !unmapped1 && aln1.ref_name == aln2.ref_name {
            "=".to_string()
        } else {
            aln2.ref_name.clone()
        };
        aln1.mate_position = aln2.position;
    }

    if unmapped1 {
        aln2.mate_ref_name = "*".to_string();
        aln2.mate_position = 0;
    } else {
        aln2.mate_ref_name = if !unmapped2 && aln1.ref_name == aln2.ref_name {
            "=".to_string()
        } else {
            aln1.ref_name.clone()
        };
        aln2.mate_position = aln1.position;
    }

    let template_length = if !unmapped1 && !unmapped2 && aln1.ref_name == aln2.ref_name {
        let start1 = aln1.position;
        let end1 = aln1.position + read1.sequence.len() as u64 - 1;
        let start2 = aln2.position;
        let end2 = aln2.position + read2.sequence.len() as u64 - 1;
        let left = start1.min(start2);
        let right = end1.max(end2);
        right as i64 - left as i64 + 1
    } else {
        0
    };

    if template_length != 0 {
        if aln1.position <= aln2.position {
            aln1.template_length = template_length;
            aln2.template_length = -template_length;
        } else {
            aln1.template_length = -template_length;
            aln2.template_length = template_length;
        }
    } else {
        aln1.template_length = 0;
        aln2.template_length = 0;
    }

    if pairing.is_proper_pair(aln1, aln2, template_length) {
        aln1.flag |= 0x2;
        aln2.flag |= 0x2;
    } else {
        aln1.flag &= !0x2;
        aln2.flag &= !0x2;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    use tempfile::NamedTempFile;

    fn read_with_len(name: &str, len: usize) -> Read {
        Read {
            name: name.to_string(),
            sequence: vec![b'A'; len],
            quality: vec![b'I'; len],
        }
    }

    fn mapped_alignment(read: &Read, ref_name: &str, position: u64, flag: u16) -> Alignment {
        Alignment {
            read_name: read.name.clone(),
            flag,
            ref_name: ref_name.to_string(),
            position,
            mapq: 60,
            cigar: format!("{}M", read.sequence.len()),
            mate_ref_name: "*".to_string(),
            mate_position: 0,
            template_length: 0,
            sequence: read.sequence.clone(),
            quality: read.quality.clone(),
            score: 100,
        }
    }

    fn write_fastq(file: &mut NamedTempFile, records: &[(&str, &str, &str)]) {
        for (name, seq, qual) in records {
            writeln!(file, "@{name}").unwrap();
            writeln!(file, "{seq}").unwrap();
            writeln!(file, "+").unwrap();
            writeln!(file, "{qual}").unwrap();
        }
        file.flush().unwrap();
    }

    #[test]
    fn sets_proper_pair_for_fr_within_insert_range() {
        let read1 = read_with_len("read/1", 100);
        let read2 = read_with_len("read/2", 100);
        let mut aln1 = mapped_alignment(&read1, "chr1", 100, 0);
        let mut aln2 = mapped_alignment(&read2, "chr1", 250, 0x10);
        let pairing = PairingHeuristic {
            min_insert_size: 0,
            max_insert_size: 500,
        };

        apply_pair_metadata(&mut aln1, &mut aln2, &read1, &read2, "read", &pairing);

        assert!(aln1.flag & 0x2 != 0);
        assert!(aln2.flag & 0x2 != 0);
    }

    #[test]
    fn does_not_set_proper_pair_for_wrong_orientation() {
        let read1 = read_with_len("read/1", 100);
        let read2 = read_with_len("read/2", 100);
        let mut aln1 = mapped_alignment(&read1, "chr1", 100, 0x10);
        let mut aln2 = mapped_alignment(&read2, "chr1", 250, 0);
        let pairing = PairingHeuristic {
            min_insert_size: 0,
            max_insert_size: 500,
        };

        apply_pair_metadata(&mut aln1, &mut aln2, &read1, &read2, "read", &pairing);

        assert!(aln1.flag & 0x2 == 0);
        assert!(aln2.flag & 0x2 == 0);
    }

    #[test]
    fn does_not_set_proper_pair_for_insert_size_out_of_range() {
        let read1 = read_with_len("read/1", 100);
        let read2 = read_with_len("read/2", 100);
        let mut aln1 = mapped_alignment(&read1, "chr1", 100, 0);
        let mut aln2 = mapped_alignment(&read2, "chr1", 2000, 0x10);
        let pairing = PairingHeuristic {
            min_insert_size: 0,
            max_insert_size: 500,
        };

        apply_pair_metadata(&mut aln1, &mut aln2, &read1, &read2, "read", &pairing);

        assert!(aln1.flag & 0x2 == 0);
        assert!(aln2.flag & 0x2 == 0);
    }

    #[test]
    fn does_not_set_proper_pair_when_mate_unmapped() {
        let read1 = read_with_len("read/1", 100);
        let read2 = read_with_len("read/2", 100);
        let mut aln1 = mapped_alignment(&read1, "chr1", 100, 0);
        let mut aln2 = Alignment::unmapped(&read2);
        let pairing = PairingHeuristic {
            min_insert_size: 0,
            max_insert_size: 500,
        };

        apply_pair_metadata(&mut aln1, &mut aln2, &read1, &read2, "read", &pairing);

        assert!(aln1.flag & 0x2 == 0);
        assert!(aln2.flag & 0x2 == 0);
    }

    #[test]
    fn read_chunk_splits_reads() {
        let mut fastq = NamedTempFile::new().unwrap();
        write_fastq(
            &mut fastq,
            &[
                ("r1", "ACGT", "!!!!"),
                ("r2", "TGCA", "####"),
                ("r3", "GATT", "$$$$"),
            ],
        );

        let mut reader = FastqReader::open(fastq.path()).unwrap();
        let first = read_chunk(&mut reader, 2).unwrap();
        assert_eq!(first.len(), 2);
        assert_eq!(first[0].name, "r1");
        assert_eq!(first[1].name, "r2");

        let second = read_chunk(&mut reader, 2).unwrap();
        assert_eq!(second.len(), 1);
        assert_eq!(second[0].name, "r3");

        let third = read_chunk(&mut reader, 2).unwrap();
        assert!(third.is_empty());
    }

    #[test]
    fn read_chunk_pairs_preserves_pair_names() {
        let mut fastq1 = NamedTempFile::new().unwrap();
        let mut fastq2 = NamedTempFile::new().unwrap();

        write_fastq(
            &mut fastq1,
            &[
                ("pair1/1", "ACGT", "!!!!"),
                ("pair2/1", "TGCA", "####"),
                ("pair3/1", "GATT", "$$$$"),
            ],
        );
        write_fastq(
            &mut fastq2,
            &[
                ("pair1/2", "ACGT", "!!!!"),
                ("pair2/2", "TGCA", "####"),
                ("pair3/2", "GATT", "$$$$"),
            ],
        );

        let mut reader1 = FastqReader::open(fastq1.path()).unwrap();
        let mut reader2 = FastqReader::open(fastq2.path()).unwrap();

        let first = read_chunk_pairs(&mut reader1, &mut reader2, 2).unwrap();
        assert_eq!(first.len(), 2);
        assert_eq!(first[0].pair_name, "pair1");
        assert_eq!(first[1].pair_name, "pair2");

        let second = read_chunk_pairs(&mut reader1, &mut reader2, 2).unwrap();
        assert_eq!(second.len(), 1);
        assert_eq!(second[0].pair_name, "pair3");

        let third = read_chunk_pairs(&mut reader1, &mut reader2, 2).unwrap();
        assert!(third.is_empty());
    }
}
