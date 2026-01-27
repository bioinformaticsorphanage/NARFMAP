//! Input/Output handling for genomic data formats

pub mod fastq;
pub mod sam;

pub use fastq::FastqReader;
pub use sam::SamWriter;
