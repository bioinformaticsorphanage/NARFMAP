use std::collections::HashMap;
use std::path::Path;
use anyhow::{Result, anyhow};
use log::{debug, info, warn};

use crate::reference::liftover::LiftCode;

/// Extend table for high-frequency k-mer positions
/// This stores actual genomic positions for k-mers that exceed frequency thresholds
#[derive(Debug)]
pub struct ExtendTable {
    /// Raw extend table data
    data: Vec<ExtendTableRecord>,
    /// Index by interval ID to table ranges (start_index, length)
    interval_index: HashMap<u64, (usize, usize)>,
    /// Configuration for extend table behavior
    config: ExtendTableConfig,
    /// Statistics for monitoring extend table usage
    stats: ExtendTableStats,
}

/// Individual extend table record following DRAGMAP's format
#[derive(Debug, Clone, Copy)]
pub struct ExtendTableRecord {
    value: u64,
}

/// Configuration for extend table behavior
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ExtendTableConfig {
    /// Enable extend table support (default: true)
    pub enabled: bool,
    /// Minimum frequency to use extend table vs. HIFREQ (default: 256)
    pub min_frequency_to_extend: u32,
    /// Maximum positions per extend table interval (default: 10000)
    pub max_interval_positions: u32,
    /// Extend table memory limit in MB (default: 512)
    pub memory_limit_mb: usize,
    /// Whether to enable interval compression
    pub enable_compression: bool,
}

/// Statistics for extend table usage
#[derive(Debug, Default, Clone)]
pub struct ExtendTableStats {
    pub total_intervals: u64,
    pub total_positions: u64,
    pub compressed_intervals: u64,
    pub memory_usage_bytes: u64,
    pub cache_hits: u64,
    pub cache_misses: u64,
}

impl Default for ExtendTableConfig {
    fn default() -> Self {
        Self {
            enabled: true,
            min_frequency_to_extend: 256,
            max_interval_positions: 10000,
            memory_limit_mb: 512,
            enable_compression: true,
        }
    }
}

impl ExtendTableRecord {
    /// Create a new extend table record
    pub fn new(position: u32, is_reverse_complement: bool, lift_code: LiftCode, lift_group: u32) -> Self {
        let mut value = 0u64;
        value |= position as u64;                                    // Position[31:0]
        value |= (is_reverse_complement as u64) << 32;              // RC[32]
        value |= ((lift_code as u64) & 0x3) << 33;                  // LiftCode[34:33]
        value |= ((lift_group as u64) & 0x1FFFFFFF) << 35;          // LiftGroup[63:35]
        Self { value }
    }

    /// Get genomic position
    pub fn position(&self) -> u32 {
        (self.value & 0xFFFFFFFF) as u32
    }

    /// Check if this is a reverse complement match
    pub fn is_reverse_complement(&self) -> bool {
        (self.value >> 32) & 1 == 1
    }

    /// Get liftover code
    pub fn lift_code(&self) -> LiftCode {
        LiftCode::from(((self.value >> 33) & 0x3) as u8)
    }

    /// Get liftover group ID
    pub fn lift_group(&self) -> u32 {
        ((self.value >> 35) & 0x1FFFFFFF) as u32
    }

    /// Get raw 64-bit value
    pub fn raw(&self) -> u64 {
        self.value
    }
}

impl ExtendTable {
    /// Create a new extend table with the given configuration
    pub fn new(config: ExtendTableConfig) -> Self {
        Self {
            data: Vec::new(),
            interval_index: HashMap::new(),
            config,
            stats: ExtendTableStats::default(),
        }
    }

    /// Add an interval of positions to the extend table
    pub fn add_interval(&mut self, interval_id: u64, positions: Vec<ExtendTableRecord>) -> Result<()> {
        if positions.len() > self.config.max_interval_positions as usize {
            warn!("Interval {} has {} positions, exceeding limit of {}. Truncating.",
                  interval_id, positions.len(), self.config.max_interval_positions);
        }

        let start_index = self.data.len();
        let actual_length = positions.len().min(self.config.max_interval_positions as usize);
        
        // Add positions to the main data vector
        self.data.extend_from_slice(&positions[..actual_length]);
        
        // Update index
        self.interval_index.insert(interval_id, (start_index, actual_length));
        
        // Update statistics
        self.stats.total_intervals += 1;
        self.stats.total_positions += actual_length as u64;
        self.stats.memory_usage_bytes += (actual_length * std::mem::size_of::<ExtendTableRecord>()) as u64;
        
        debug!("Added interval {} with {} positions starting at index {}", 
               interval_id, actual_length, start_index);
        
        Ok(())
    }

    /// Query positions for a given interval ID
    pub fn query_interval(&mut self, interval_id: u64) -> Option<&[ExtendTableRecord]> {
        if let Some(&(start_index, length)) = self.interval_index.get(&interval_id) {
            self.stats.cache_hits += 1;
            Some(&self.data[start_index..start_index + length])
        } else {
            self.stats.cache_misses += 1;
            None
        }
    }

    /// Get extend table statistics
    pub fn get_stats(&self) -> &ExtendTableStats {
        &self.stats
    }

    /// Check if extend table is enabled and within memory limits
    pub fn is_usable(&self) -> bool {
        self.config.enabled && 
        (self.stats.memory_usage_bytes / (1024 * 1024)) <= self.config.memory_limit_mb as u64
    }

    /// Save extend table to disk
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        use std::fs::File;
        use std::io::{BufWriter, Write};
        
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        
        // Write header with metadata
        let header = ExtendTableHeader {
            version: 1,
            config: self.config.clone(),
            stats: self.stats.clone(),
            data_length: self.data.len() as u64,
            index_length: self.interval_index.len() as u64,
        };
        
        // Serialize header as JSON for now (could use binary format later)
        let header_json = serde_json::to_string(&header)?;
        let header_len = header_json.len() as u32;
        
        writer.write_all(&header_len.to_le_bytes())?;
        writer.write_all(header_json.as_bytes())?;
        
        // Write data records
        for record in &self.data {
            writer.write_all(&record.value.to_le_bytes())?;
        }
        
        // Write index
        let index_len = self.interval_index.len() as u32;
        writer.write_all(&index_len.to_le_bytes())?;
        
        for (&interval_id, &(start_index, length)) in &self.interval_index {
            writer.write_all(&interval_id.to_le_bytes())?;
            writer.write_all(&(start_index as u64).to_le_bytes())?;
            writer.write_all(&(length as u64).to_le_bytes())?;
        }
        
        writer.flush()?;
        info!("Saved extend table with {} intervals and {} positions", 
              self.stats.total_intervals, self.stats.total_positions);
        
        Ok(())
    }

    /// Load extend table from disk
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        use std::fs::File;
        use std::io::{BufReader, Read};
        
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);
        
        // Read header length
        let mut header_len_bytes = [0u8; 4];
        reader.read_exact(&mut header_len_bytes)?;
        let header_len = u32::from_le_bytes(header_len_bytes) as usize;
        
        // Read header
        let mut header_bytes = vec![0u8; header_len];
        reader.read_exact(&mut header_bytes)?;
        let header_json = String::from_utf8(header_bytes)?;
        let header: ExtendTableHeader = serde_json::from_str(&header_json)?;
        
        // Read data records
        let mut data = Vec::with_capacity(header.data_length as usize);
        for _ in 0..header.data_length {
            let mut record_bytes = [0u8; 8];
            reader.read_exact(&mut record_bytes)?;
            let value = u64::from_le_bytes(record_bytes);
            data.push(ExtendTableRecord { value });
        }
        
        // Read index length
        let mut index_len_bytes = [0u8; 4];
        reader.read_exact(&mut index_len_bytes)?;
        let index_len = u32::from_le_bytes(index_len_bytes) as usize;
        
        // Read index
        let mut interval_index = HashMap::with_capacity(index_len);
        for _ in 0..index_len {
            let mut interval_id_bytes = [0u8; 8];
            let mut start_index_bytes = [0u8; 8];
            let mut length_bytes = [0u8; 8];
            
            reader.read_exact(&mut interval_id_bytes)?;
            reader.read_exact(&mut start_index_bytes)?;
            reader.read_exact(&mut length_bytes)?;
            
            let interval_id = u64::from_le_bytes(interval_id_bytes);
            let start_index = u64::from_le_bytes(start_index_bytes) as usize;
            let length = u64::from_le_bytes(length_bytes) as usize;
            
            interval_index.insert(interval_id, (start_index, length));
        }
        
        info!("Loaded extend table with {} intervals and {} positions", 
              header.stats.total_intervals, header.stats.total_positions);
        
        Ok(Self {
            data,
            interval_index,
            config: header.config,
            stats: header.stats,
        })
    }

    /// Clear all data and reset statistics
    pub fn clear(&mut self) {
        self.data.clear();
        self.interval_index.clear();
        self.stats = ExtendTableStats::default();
    }

    /// Get total memory usage in bytes
    pub fn memory_usage(&self) -> u64 {
        let data_size = self.data.len() * std::mem::size_of::<ExtendTableRecord>();
        let index_size = self.interval_index.len() * (std::mem::size_of::<u64>() + std::mem::size_of::<(usize, usize)>());
        (data_size + index_size) as u64
    }
}

/// Header for extend table serialization
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
struct ExtendTableHeader {
    version: u32,
    config: ExtendTableConfig,
    stats: ExtendTableStats,
    data_length: u64,
    index_length: u64,
}

// Make ExtendTableStats serializable
impl serde::Serialize for ExtendTableStats {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        use serde::ser::SerializeStruct;
        let mut state = serializer.serialize_struct("ExtendTableStats", 6)?;
        state.serialize_field("total_intervals", &self.total_intervals)?;
        state.serialize_field("total_positions", &self.total_positions)?;
        state.serialize_field("compressed_intervals", &self.compressed_intervals)?;
        state.serialize_field("memory_usage_bytes", &self.memory_usage_bytes)?;
        state.serialize_field("cache_hits", &self.cache_hits)?;
        state.serialize_field("cache_misses", &self.cache_misses)?;
        state.end()
    }
}

impl<'de> serde::Deserialize<'de> for ExtendTableStats {
    fn deserialize<D>(deserializer: D) -> Result<ExtendTableStats, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        use serde::de::{self, Deserialize, Deserializer, MapAccess, Visitor};
        use std::fmt;

        #[derive(serde::Deserialize)]
        #[serde(field_identifier, rename_all = "snake_case")]
        enum Field {
            TotalIntervals,
            TotalPositions,
            CompressedIntervals,
            MemoryUsageBytes,
            CacheHits,
            CacheMisses,
        }

        struct ExtendTableStatsVisitor;

        impl<'de> Visitor<'de> for ExtendTableStatsVisitor {
            type Value = ExtendTableStats;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("struct ExtendTableStats")
            }

            fn visit_map<V>(self, mut map: V) -> Result<ExtendTableStats, V::Error>
            where
                V: MapAccess<'de>,
            {
                let mut total_intervals = None;
                let mut total_positions = None;
                let mut compressed_intervals = None;
                let mut memory_usage_bytes = None;
                let mut cache_hits = None;
                let mut cache_misses = None;

                while let Some(key) = map.next_key()? {
                    match key {
                        Field::TotalIntervals => {
                            if total_intervals.is_some() {
                                return Err(de::Error::duplicate_field("total_intervals"));
                            }
                            total_intervals = Some(map.next_value()?);
                        }
                        Field::TotalPositions => {
                            if total_positions.is_some() {
                                return Err(de::Error::duplicate_field("total_positions"));
                            }
                            total_positions = Some(map.next_value()?);
                        }
                        Field::CompressedIntervals => {
                            if compressed_intervals.is_some() {
                                return Err(de::Error::duplicate_field("compressed_intervals"));
                            }
                            compressed_intervals = Some(map.next_value()?);
                        }
                        Field::MemoryUsageBytes => {
                            if memory_usage_bytes.is_some() {
                                return Err(de::Error::duplicate_field("memory_usage_bytes"));
                            }
                            memory_usage_bytes = Some(map.next_value()?);
                        }
                        Field::CacheHits => {
                            if cache_hits.is_some() {
                                return Err(de::Error::duplicate_field("cache_hits"));
                            }
                            cache_hits = Some(map.next_value()?);
                        }
                        Field::CacheMisses => {
                            if cache_misses.is_some() {
                                return Err(de::Error::duplicate_field("cache_misses"));
                            }
                            cache_misses = Some(map.next_value()?);
                        }
                    }
                }

                Ok(ExtendTableStats {
                    total_intervals: total_intervals.unwrap_or(0),
                    total_positions: total_positions.unwrap_or(0),
                    compressed_intervals: compressed_intervals.unwrap_or(0),
                    memory_usage_bytes: memory_usage_bytes.unwrap_or(0),
                    cache_hits: cache_hits.unwrap_or(0),
                    cache_misses: cache_misses.unwrap_or(0),
                })
            }
        }

        deserializer.deserialize_struct("ExtendTableStats", &["total_intervals", "total_positions", "compressed_intervals", "memory_usage_bytes", "cache_hits", "cache_misses"], ExtendTableStatsVisitor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn test_extend_table_record() {
        let record = ExtendTableRecord::new(
            12345,                    // position
            true,                     // reverse complement
            LiftCode::None,           // lift code
            67890                     // lift group
        );

        assert_eq!(record.position(), 12345);
        assert!(record.is_reverse_complement());
        assert_eq!(record.lift_code(), LiftCode::None);
        assert_eq!(record.lift_group(), 67890);
    }

    #[test]
    fn test_extend_table_operations() {
        let mut extend_table = ExtendTable::new(ExtendTableConfig::default());

        // Add an interval
        let positions = vec![
            ExtendTableRecord::new(100, false, LiftCode::None, 1),
            ExtendTableRecord::new(200, true, LiftCode::None, 1),
            ExtendTableRecord::new(300, false, LiftCode::Alt, 2),
        ];

        extend_table.add_interval(1, positions.clone()).unwrap();

        // Query the interval
        let result = extend_table.query_interval(1).unwrap();
        assert_eq!(result.len(), 3);
        assert_eq!(result[0].position(), 100);
        assert_eq!(result[1].position(), 200);
        assert_eq!(result[2].position(), 300);

        // Query non-existent interval
        assert!(extend_table.query_interval(999).is_none());

        // Check statistics
        let stats = extend_table.get_stats();
        assert_eq!(stats.total_intervals, 1);
        assert_eq!(stats.total_positions, 3);
        assert_eq!(stats.cache_hits, 1);
        assert_eq!(stats.cache_misses, 1);
    }

    #[test]
    fn test_extend_table_serialization() {
        let temp_dir = TempDir::new().unwrap();
        let file_path = temp_dir.path().join("test_extend_table.bin");

        // Create and populate extend table
        let mut extend_table = ExtendTable::new(ExtendTableConfig::default());
        let positions = vec![
            ExtendTableRecord::new(100, false, LiftCode::None, 1),
            ExtendTableRecord::new(200, true, LiftCode::Alt, 2),
        ];
        extend_table.add_interval(42, positions).unwrap();

        // Save to disk
        extend_table.save(&file_path).unwrap();

        // Load from disk
        let loaded_table = ExtendTable::load(&file_path).unwrap();

        // Verify data integrity
        let result = loaded_table.interval_index.get(&42).unwrap();
        assert_eq!(result.1, 2); // length

        assert_eq!(loaded_table.data.len(), 2);
        assert_eq!(loaded_table.data[0].position(), 100);
        assert_eq!(loaded_table.data[1].position(), 200);
    }
}