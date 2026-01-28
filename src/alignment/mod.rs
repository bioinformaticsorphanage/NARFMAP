//! Alignment algorithms

pub mod aligner;
pub mod mapper;
pub mod mapq;

pub use aligner::Aligner;
pub use mapper::SeedMapper;
pub use mapq::{compute_mapq, MAPQ_MAX};
