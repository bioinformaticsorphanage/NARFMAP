# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NARFMAP is a Rust rewrite of the DRAGMAP aligner, which is a high-performance DNA sequence aligner. The project is currently in active transition from a C++ codebase to pure Rust implementation.

**Current Status**: Modern CLI interface completed, hash table generation implemented, Smith-Waterman alignment implemented, alt-aware mapping with liftover support implemented. Approximately 85-90% DRAGMAP feature parity achieved.

## Build Commands

### Rust Development (Primary)
```bash
# Build the project
cargo build
cargo build --release

# Run tests
cargo test

# Run the modern CLI
cargo run -- --help

# Index a reference genome
cargo run -- index reference.fa

# Align reads (single-end)
cargo run -- map reference.fa reads.fastq

# Align reads (paired-end)
cargo run -- map reference.fa reads_1.fastq reads_2.fastq

# Show file statistics
cargo run -- stats reads.fastq

# Validate file integrity
cargo run -- check reads.fastq

# Inspect an existing index
cargo run -- inspect reference.fa.index

# Generate configuration file
cargo run -- config --output my_config.toml

# Run with logging
RUST_LOG=debug cargo run -- <command>
```

### Legacy C++ Build (Fallback)
```bash
# Build C++ version (for testing/comparison)
make
make install

# Build without tests
HAS_GTEST=0 make

# Clean build
make clean
```

### Nix (Reproducible Builds)
```bash
nix build
```

## Architecture

### Command Structure
- **index**: Build reference genome index with k-mer hash tables
- **map**: Align FASTQ reads to reference genome (single-end or paired-end)
- **stats**: Display comprehensive statistics about sequence files
- **check**: Validate file integrity and format correctness
- **config**: Generate configuration files with presets for different use cases
- **inspect**: Examine existing index files and show detailed information

### Rust Module Organization
```
src/
├── main.rs          # CLI entry point and command dispatch
├── cli.rs           # Modern CLI interface with clap, colored output, and progress bars
├── align/           # Alignment algorithms with Smith-Waterman implementation
│   ├── mod.rs       # Main aligner with seed-and-extend pipeline
│   ├── smith_waterman.rs  # Full Smith-Waterman with DP matrix
│   ├── stats.rs     # Alignment statistics tracking
│   ├── sam.rs       # SAM format output generation
│   └── paired.rs    # Paired-end alignment and mate rescue
├── config/          # Configuration management with DRAGMAP compatibility
├── hashtable/       # Hash table implementation with k-mer indexing
├── io/              # FASTA/FASTQ file handling
├── reference/       # Reference genome processing with liftover support
│   └── liftover.rs  # Alt-aware mapping and liftover functionality
└── utils/           # Utility functions
```

### Key Dependencies
- **CLI**: `clap` with derive features, `colored` for terminal output, `indicatif` for progress bars
- **Bioinformatics**: `seq_io`, `bitnuc`, `noodles`
- **Concurrency**: `rayon`, `crossbeam`, `dashmap`
- **I/O**: `flate2`, `memmap2`
- **Error Handling**: `anyhow`, `thiserror`
- **Configuration**: `toml` for config files, `serde` for serialization
- **Time/Date**: `chrono` for SAM headers and statistics
- **System**: `atty` for terminal detection, `num_cpus` for thread count

## Development Notes

### Current Implementation Status
- ✅ Modern CLI interface with intuitive commands and colored output
- ✅ Basic I/O modules for FASTA/FASTQ with progress bars
- ✅ Logging and error handling setup
- ✅ Configuration system with TOML files and presets
- ✅ CRC-based k-mer hashing with 2-bit nucleotide encoding
- ✅ Reference sequence processing (FASTA loading, k-mer extraction)
- ✅ Hash table generation with bucket-based storage and serialization
- ✅ Smith-Waterman alignment with dynamic programming matrix and banded optimization
- ✅ Alt-aware mapping with liftover support for alternative contigs
- ✅ Liftover group management for primary/alt contig relationships
- ✅ Mate rescue scanning for paired-end reads within insert size intervals
- ✅ Dynamic seed extension based on k-mer frequency thresholds
- ✅ Seed-and-extend alignment pipeline with clustering
- ✅ Liftover-based MAPQ calculation with group scoring
- ✅ Empirical insert size distribution modeling and updates
- ✅ Comprehensive alignment statistics tracking
- ✅ Unit and integration test framework

### Remaining High-Priority DRAGMAP Features
1. **Reference masking**: Support for masked regions in alt contigs
2. **Extension chaining**: Sophisticated chaining logic for seed clusters
3. **Split alignment discovery**: Detection and scoring for large indels
4. **Hash table probing**: Linear probing and chaining for collision handling
5. **Extend table support**: High-frequency k-mer handling

### Medium-Priority Features
1. **Wavefront steering**: Smith-Waterman optimization
2. **IUB codes injection**: Multi-nucleotide codes for isolated SNVs
3. **Frequency-based seed extension filtering**: Advanced filtering logic
4. **Supplementary alignment generation**: Beyond SAM flag support
5. **Global alignment option**: End-to-end read alignment

### Legacy C++ Reference
The `src/include/` and `src/lib/` directories contain the original C++ implementation that can serve as a reference for the Rust rewrite. Key components:
- Hash table algorithms in `src/lib/hashtable/`
- Alignment logic in `src/lib/align/`
- Reference processing in `src/lib/reference/`

### Testing Strategy
- Legacy C++ tests use GoogleTest framework
- Rust tests should use standard `#[test]` functions
- Use `criterion` for benchmarking performance-critical code
- Test data available in `data/` directory

### Development Workflow
1. Implement core algorithms one module at a time
2. Write tests alongside implementation
3. Use `cargo test` for validation
4. Profile with `cargo bench` for performance
5. Compare output with legacy C++ implementation for correctness
6. **Make small, focused commits**: Commit related changes together in logical groups for easier review and debugging

### Performance Considerations
- Use `rayon` for parallel processing
- Memory-map large files with `memmap2`
- Optimize hot paths identified through profiling
- Target performance parity with C++ implementation

### DRAGMAP Compatibility Status
NARFMAP now implements approximately **90-92% of DRAGMAP's core features**:

**✅ Implemented (Core Features)**
- k-mer hashing with CRC polynomials (100% compatible)
- Hash table generation with bucket storage
- Hash table probing and chaining for collision handling (block-constrained)
- Extend table support for high-frequency k-mers (≥256 occurrences)
- Smith-Waterman alignment with banded optimization
- Dynamic seed extension (matches DRAGMAP's MAX_EXTENSION_BASES)
- Seed clustering and alignment candidate evaluation
- MAPQ calculation with primary/secondary score comparison
- DRAGMAP-compatible scoring matrices and penalties
- Alt-aware mapping and liftover groups (major DRAGMAP differentiator)
- Mate rescue scanning for paired-end reads
- Liftover-based MAPQ calculation with group scoring
- Empirical insert size distribution modeling and updates
- Paired-end orientation detection and validation

**🔄 Partially Implemented**
- SAM format output (comprehensive implementation, missing some advanced features)
- Reference processing (FASTA loading complete, needs masking support)

**❌ Not Yet Implemented**
- Split alignment discovery for structural variants
- Wavefront steering optimizations
- Reference masking for alt contigs
- Extension chaining logic for seed clusters
- IUB codes injection for isolated SNVs

### Key Implementation Insights
- DRAGMAP uses complex template-based C++ with sophisticated type hierarchies
- NARFMAP simplifies this with Rust's type system while maintaining algorithmic fidelity  
- Smith-Waterman implementation follows classic DP approach with DRAGMAP scoring
- Dynamic seed extension successfully mimics DRAGMAP's frequency-based thresholds
- Current implementation prioritizes correctness over performance optimizations