# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

NARFMAP is a Rust rewrite of the DRAGMAP aligner, which is a high-performance DNA sequence aligner. The project is currently in active transition from a C++ codebase to pure Rust implementation.

**Current Status**: CLI interface completed, core algorithms (hash table generation, alignment) are stubs that need implementation.

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
├── align/           # Alignment algorithms (TODO: implement)
├── config/          # Configuration management
├── hashtable/       # Hash table implementation (TODO: implement)
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

## Development Notes

### Current Implementation Status
- ✅ CLI interface and argument parsing
- ✅ Basic I/O modules for FASTA/FASTQ
- ✅ Logging and error handling setup
- ✅ Configuration system with DRAGMAP-compatible parameters
- ✅ CRC-based k-mer hashing with 2-bit nucleotide encoding
- ✅ Reference sequence processing (FASTA loading, k-mer extraction)
- 🔄 Hash table generation (core structures implemented, serialization pending)
- ❌ Alignment algorithms (returns `unimplemented!()`)
- ✅ Unit and integration test framework

### Critical Implementation Needs
1. **Hash table generation**: Core k-mer indexing functionality
2. **Alignment algorithms**: Seed-and-extend alignment logic
3. **Reference processing**: Efficient reference genome handling
4. **Test coverage**: Unit and integration tests for Rust modules

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

### Project Status Notes
- The remaining tasks (linear probing, reference names, alignment extension) are optimizations and enhancements rather than core functionality fixes.