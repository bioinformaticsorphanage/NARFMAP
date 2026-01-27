# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Add prek pre-commit hooks (fmt, clippy, check, test) with fail-fast for agent workflows

## [1.4.2] - 2024-07-07

### Fixed

- Add a missing import to fix Bioconda build [19e3063](https://github.com/bioinformaticsorphanage/NARFMAP/commit/19e3063)

## [1.4.1] - 2024-07-05

### Changed

- Refactor flake to use flake-parts
- Add treefmt

### Fixed

- Allow to statically compile dragmap by @gambalab in https://github.com/bioinformaticsorphanage/NARFMAP/pull/12

## [1.4.0] - 2023-12-07

### Added

- CI Action to build (#1)
- Add v2 FASTA masks from Illumina from [ckandoth](https://github.com/ckandoth). ([DRAGMAP #32](https://github.com/Illumina/DRAGMAP/pull/32))
- Add nix flake (#4)
- Add initial Rust infrastructure
- Add new Bioconda recipe (#5)

### Fixed

- Remove useless space in @PG from [lindenb](https://github.com/lindenb). ([DRAGMAP #14](https://github.com/Illumina/DRAGMAP/pull/14))
- Fix overflow TableCount from [tomofumisaka](https://github.com/tomofumisaka). ([DRAGMAP #55](https://github.com/Illumina/DRAGMAP/pull/55))
- Fix overflow number of reference contig from [tomofumisaka](https://github.com/tomofumisaka). ([DRAGMAP #56](https://github.com/Illumina/DRAGMAP/pull/56))
- Fix "Assertion `inputR2.eof()' failed." error from [lgruen](https://github.com/lgruen). ([DRAGMAP #53](https://github.com/Illumina/DRAGMAP/pull/53))
- Fix intel compiler build break from [fo40225](https://github.com/fo40225). ([DRAGMAP #48](https://github.com/Illumina/DRAGMAP/pull/48))

### Changed

- Use gnu99 instead of c99from [lgruen](https://github.com/lgruen). ([DRAGMAP #10](https://github.com/Illumina/DRAGMAP/pull/10))
- Add treefmt (#6)

## [1.3.1] - 2022-12-08

### Fixed

- mapq0 fix

[unreleased]: https://github.com/Emiller88/NARFMAP/compare/v1.4.2...HEAD
[1.4.2]: https://github.com/bioinformaticsorphanage/NARFMAP/compare/1.4.1...1.4.2
[1.4.1]: https://github.com/bioinformaticsorphanage/NARFMAP/compare/1.4.0...1.4.1
[1.4.0]: https://github.com/Emiller88/NARFMAP/compare/v1.3.1...v1.4.0
[1.3.1]: https://github.com/Emiller88/NARFMAP/compare/v1.3.0...v1.3.1
