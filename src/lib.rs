// Library interface for NARFMAP
// This exposes the core functionality for integration testing and potential library use

pub mod align;
pub mod config;
pub mod hashtable;
pub mod io;
pub mod reference;
pub mod utils;

// Re-export commonly used types for easier access
pub use align::{Aligner, AlignmentResult};
pub use config::{Config, AlignmentConfig, HashTableConfig};
pub use hashtable::{HashTableBuilder, HashTableQuery};
pub use io::{fastq::FastqReader, sequence::Sequence};