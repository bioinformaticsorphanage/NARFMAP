# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NARFMAP is a Rust rewrite of the DRAGMAP aligner, which is a high-performance DNA sequence aligner. The project is currently in active transition from a C++ codebase to pure Rust implementation.

**Current Status**: CLI interface completed, hash table generation implemented, Smith-Waterman alignment implemented. Approximately 75-80% DRAGMAP feature parity achieved.

## Build Commands

### Rust Development (Primary)
```bash
# Build the project
cargo build
cargo build --release

# Run tests
cargo test

# Run the CLI
cargo run -- --help
cargo run -- build-hash-table -r reference.fa
cargo run -- align -r ref_dir -1 reads.fastq
cargo run -- info -i sequences.fa

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
- **build-hash-table**: Generate reference hash tables for alignment
- **align**: Align FASTQ reads against reference hash table
- **info**: Display information about FASTA/FASTQ files

### Rust Module Organization
```
src/
├── main.rs          # CLI entry point with clap
├── align/           # Alignment algorithms with Smith-Waterman implementation
│   ├── mod.rs       # Main aligner with seed-and-extend pipeline
│   ├── smith_waterman.rs  # Full Smith-Waterman with DP matrix
│   ├── stats.rs     # Alignment statistics tracking
│   └── sam.rs       # SAM format output generation
├── config/          # Configuration management with DRAGMAP compatibility
├── hashtable/       # Hash table implementation with k-mer indexing
├── io/              # FASTA/FASTQ file handling
├── reference/       # Reference genome processing
└── utils/           # Utility functions
```

### Key Dependencies
- **CLI**: `clap` with derive features
- **Bioinformatics**: `seq_io`, `bitnuc`, `noodles`
- **Concurrency**: `rayon`, `crossbeam`, `dashmap`
- **I/O**: `flate2`, `memmap2`
- **Error Handling**: `anyhow`, `thiserror`
- **Time/Date**: `chrono` for SAM headers and statistics

## Development Notes

### Current Implementation Status
- ✅ CLI interface and argument parsing
- ✅ Basic I/O modules for FASTA/FASTQ
- ✅ Logging and error handling setup
- ✅ Configuration system with DRAGMAP-compatible parameters
- ✅ CRC-based k-mer hashing with 2-bit nucleotide encoding
- ✅ Reference sequence processing (FASTA loading, k-mer extraction)
- ✅ Hash table generation with bucket-based storage and serialization
- ✅ Smith-Waterman alignment with dynamic programming matrix
- ✅ Dynamic seed extension based on k-mer frequency thresholds
- ✅ Seed-and-extend alignment pipeline with clustering
- ✅ MAPQ calculation with primary/secondary score comparison
- ✅ Comprehensive alignment statistics tracking
- ✅ Unit and integration test framework

### Remaining High-Priority DRAGMAP Features
1. **Alt-aware mapping**: Liftover support for alternative contigs
2. **Mate rescue scanning**: Paired-end mate recovery within insert size intervals
3. **Reference masking**: Support for masked regions in alt contigs
4. **Extension chaining**: Sophisticated chaining logic for seed clusters
5. **Insert size distribution**: Empirical modeling and updates for paired-end reads

### Medium-Priority Features
1. **Split alignment discovery**: Detection and scoring for large indels
2. **Hash table probing**: Linear probing and chaining for collision handling
3. **Extend table support**: High-frequency k-mer handling
4. **Wavefront steering**: Smith-Waterman optimization
5. **IUB codes injection**: Multi-nucleotide codes for isolated SNVs

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
NARFMAP now implements approximately **75-80% of DRAGMAP's core features**:

**✅ Implemented (Core Features)**
- k-mer hashing with CRC polynomials (100% compatible)
- Hash table generation with bucket storage
- Smith-Waterman alignment with banded optimization
- Dynamic seed extension (matches DRAGMAP's MAX_EXTENSION_BASES)
- Seed clustering and alignment candidate evaluation
- MAPQ calculation with primary/secondary score comparison
- DRAGMAP-compatible scoring matrices and penalties

**🔄 Partially Implemented**
- SAM format output (basic implementation, missing advanced features)
- Paired-end support (framework exists, needs insert size modeling)
- Reference processing (basic FASTA loading, needs masking support)

**❌ Not Yet Implemented**
- Alt-aware mapping and liftover groups (major DRAGMAP differentiator)
- Mate rescue scanning for paired-end reads
- Split alignment discovery for structural variants
- Hash table extend tables for high-frequency k-mers
- Wavefront steering optimizations

### Key Implementation Insights
- DRAGMAP uses complex template-based C++ with sophisticated type hierarchies
- NARFMAP simplifies this with Rust's type system while maintaining algorithmic fidelity  
- Smith-Waterman implementation follows classic DP approach with DRAGMAP scoring
- Dynamic seed extension successfully mimics DRAGMAP's frequency-based thresholds
- Current implementation prioritizes correctness over performance optimizations