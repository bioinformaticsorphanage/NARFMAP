//! Reference genome and hash table handling

pub mod hash_table;
pub mod sequence;

pub use hash_table::{HashTable, HashTableConfig};
pub use sequence::ReferenceSequence;
