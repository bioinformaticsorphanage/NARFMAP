// use std::path::PathBuf; // Unused for now
use serde::{Deserialize, Serialize};
use anyhow::Result;

/// Hash table generation configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HashTableConfig {
    /// Initial seed length to store in hash table (default: 21)
    pub seed_len: usize,
    /// Maximum extended seed length (default: 0, auto-calculated as seed_len + 128)
    pub max_ext_seed_len: usize,
    /// Number of positions per reference seed (default: 1.0)
    pub ref_seed_interval: f64,
    /// Maximum allowed frequency for seed matches after extension (1-256, default: 16)
    pub max_seed_freq: u32,
    /// Target seed frequency for seed extension (default: 4.0)
    pub target_seed_freq: f64,
    /// Soft seed frequency cap for thinning (default: 12.0)
    pub soft_seed_freq_cap: f64,
    /// Maximum bases to extend a seed by in one step (default: 12)
    pub max_ext_incr: usize,
    /// Number of threads for hash table generation (default: 8)
    pub num_threads: usize,
    /// Memory limit for hash table + reference in GB (default: 32)
    pub mem_limit_gb: usize,
    /// CRC polynomial index for hashing primary seeds (default: 0)
    pub crc_primary: u32,
    /// CRC polynomial index for hashing extended seeds (default: 0)
    pub crc_extended: u32,
}

impl Default for HashTableConfig {
    fn default() -> Self {
        Self {
            seed_len: 21,
            max_ext_seed_len: 0, // Will be calculated as seed_len + 128
            ref_seed_interval: 1.0,
            max_seed_freq: 16,
            target_seed_freq: 4.0,
            soft_seed_freq_cap: 12.0,
            max_ext_incr: 12,
            num_threads: 8,
            mem_limit_gb: 32,
            crc_primary: 20, // 20 bits for reasonable hash table size
            crc_extended: 54, // 54 bits for extended hashing
        }
    }
}

/// Alignment scoring configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentConfig {
    /// Score for matching bases (default: 1)
    pub match_score: i32,
    /// Penalty for mismatched bases (default: -4)
    pub mismatch_score: i32,
    /// Gap opening penalty (default: 7)
    pub gap_init_penalty: i32,
    /// Gap extension penalty (default: 1)
    pub gap_extend_penalty: i32,
    /// Unclipping score (default: 5)
    pub unclip_score: i32,
    /// Minimum alignment score (default: 22)
    pub min_score: i32,
    /// Number of worker threads for alignment (default: hardware threads)
    pub num_threads: usize,
    /// Maximum secondary alignments per read (default: 0)
    pub max_secondary_aligns: u32,
    /// Score delta for secondary alignments (default: 0)
    pub secondary_score_delta: i32,
    /// Seed length for alignment (default: 21)
    pub seed_len: usize,
    /// Step size for extracting seeds from reads (default: 1)
    pub seed_step_size: usize,
    /// Maximum distance to cluster seed hits (default: 1000)
    pub cluster_distance: u32,
}

impl Default for AlignmentConfig {
    fn default() -> Self {
        Self {
            match_score: 1,
            mismatch_score: -4,
            gap_init_penalty: 7,
            gap_extend_penalty: 1,
            unclip_score: 5,
            min_score: 22,
            num_threads: std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(4),
            max_secondary_aligns: 0,
            secondary_score_delta: 0,
            seed_len: 21,
            seed_step_size: 1,
            cluster_distance: 1000,
        }
    }
}

/// Paired-end configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PairedEndConfig {
    /// Expected mean insert size (default: 0.0, auto-detect)
    pub mean_insert_size: f64,
    /// Expected standard deviation of insert size (default: 0.0, auto-detect)
    pub stddev_insert_size: f64,
    /// Expected mean read length (default: 0, auto-detect)
    pub mean_read_len: usize,
    /// Expected orientation: 0=FR, 1=RF, 2=FF (default: 0)
    pub orientation: u32,
    /// Enable auto-detection of paired-end parameters (default: true)
    pub auto_detect_params: bool,
    /// Recent pairs for stats calculation (default: 100000)
    pub stats_sample_size: usize,
}

impl Default for PairedEndConfig {
    fn default() -> Self {
        Self {
            mean_insert_size: 0.0,
            stddev_insert_size: 0.0,
            mean_read_len: 0,
            orientation: 0, // FR orientation
            auto_detect_params: true,
            stats_sample_size: 100000,
        }
    }
}

/// I/O configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IOConfig {
    /// FASTQ quality score offset (33 or 64, default: 33)
    pub fastq_offset: u8,
    /// Character for qname suffix delimiter (default: ' ')
    pub qname_suffix_delimiter: char,
    /// Interleaved paired-end reads flag (default: false)
    pub interleaved: bool,
    /// Preserve deterministic output order (default: false)
    pub preserve_order: bool,
    /// Read group ID (default: "1")
    pub read_group_id: String,
    /// Read group sample name (default: "sample")
    pub read_group_sample: String,
}

impl Default for IOConfig {
    fn default() -> Self {
        Self {
            fastq_offset: 33,
            qname_suffix_delimiter: ' ',
            interleaved: false,
            preserve_order: false,
            read_group_id: "1".to_string(),
            read_group_sample: "sample".to_string(),
        }
    }
}

/// Main configuration structure combining all settings
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Config {
    /// Hash table configuration
    pub hashtable: HashTableConfig,
    /// Alignment configuration
    pub alignment: AlignmentConfig,
    /// Paired-end configuration
    pub paired_end: PairedEndConfig,
    /// I/O configuration
    pub io: IOConfig,
    /// Verbosity level (0=warn, 1=info, 2=debug, 3=trace)
    pub verbosity: u8,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            hashtable: HashTableConfig::default(),
            alignment: AlignmentConfig::default(),
            paired_end: PairedEndConfig::default(),
            io: IOConfig::default(),
            verbosity: 1, // info level by default
        }
    }
}

impl Config {
    /// Create a new configuration with default values
    pub fn new() -> Self {
        Self::default()
    }

    /// Load configuration from a file (JSON or TOML)
    pub fn from_file<P: AsRef<std::path::Path>>(path: P) -> Result<Self> {
        let content = std::fs::read_to_string(&path)?;
        let path_ref = path.as_ref();
        
        if let Some(ext) = path_ref.extension().and_then(|s| s.to_str()) {
            match ext.to_lowercase().as_str() {
                "json" => Ok(serde_json::from_str(&content)?),
                "toml" => Ok(toml::from_str(&content)?),
                _ => {
                    // Try JSON first, then TOML
                    serde_json::from_str(&content)
                        .or_else(|_| toml::from_str(&content))
                        .map_err(|e| anyhow::anyhow!("Failed to parse config file: {}", e))
                }
            }
        } else {
            // Try JSON first, then TOML
            serde_json::from_str(&content)
                .or_else(|_| toml::from_str(&content))
                .map_err(|e| anyhow::anyhow!("Failed to parse config file: {}", e))
        }
    }

    /// Save configuration to a file
    pub fn to_file<P: AsRef<std::path::Path>>(&self, path: P) -> Result<()> {
        let path_ref = path.as_ref();
        let content = if let Some(ext) = path_ref.extension().and_then(|s| s.to_str()) {
            match ext.to_lowercase().as_str() {
                "json" => serde_json::to_string_pretty(self)?,
                "toml" => toml::to_string_pretty(self)?,
                _ => serde_json::to_string_pretty(self)?, // Default to JSON
            }
        } else {
            serde_json::to_string_pretty(self)? // Default to JSON
        };
        
        std::fs::write(path, content)?;
        Ok(())
    }

    /// Validate configuration parameters
    pub fn validate(&mut self) -> Result<()> {
        // Auto-calculate max_ext_seed_len if not set
        if self.hashtable.max_ext_seed_len == 0 {
            self.hashtable.max_ext_seed_len = self.hashtable.seed_len + 128;
        }

        // Validate seed length
        if self.hashtable.seed_len < 8 || self.hashtable.seed_len > 32 {
            return Err(anyhow::anyhow!("Seed length must be between 8 and 32"));
        }

        // Validate seed frequency
        if self.hashtable.max_seed_freq < 1 || self.hashtable.max_seed_freq > 256 {
            return Err(anyhow::anyhow!("Max seed frequency must be between 1 and 256"));
        }

        // Validate thread counts
        if self.hashtable.num_threads == 0 {
            return Err(anyhow::anyhow!("Hash table thread count must be > 0"));
        }
        if self.alignment.num_threads == 0 {
            return Err(anyhow::anyhow!("Alignment thread count must be > 0"));
        }

        // Validate FASTQ offset
        if self.io.fastq_offset != 33 && self.io.fastq_offset != 64 {
            return Err(anyhow::anyhow!("FASTQ offset must be 33 or 64"));
        }

        // Validate paired-end orientation
        if self.paired_end.orientation > 2 {
            return Err(anyhow::anyhow!("Paired-end orientation must be 0, 1, or 2"));
        }

        Ok(())
    }

    /// Apply command-line overrides
    pub fn override_from_cli(&mut self, _args: &()) {
        // TODO: Implement CLI overrides when Cli struct is accessible
        // self.verbosity = args.verbose;
        // Additional CLI overrides can be added here as needed
    }
} 