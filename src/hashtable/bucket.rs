use crate::hashtable::hash_record::HashRecord;

/// Hash table bucket containing 8 hash records (64 bytes)
/// 
/// This matches DRAGMAP's bucket structure:
/// - 8 records per bucket (BUCKET_RECORDS_LOG2 = 3)
/// - 8 bytes per record (HASH_RECORD_BYTES_LOG2 = 3)
/// - 64 bytes per bucket (HASH_BUCKET_BYTES_LOG2 = 6)
#[derive(Clone, Copy, Debug)]
pub struct Bucket {
    records: [HashRecord; 8],
}

/// Hash table traits matching DRAGMAP constants
pub struct HashtableTraits;

impl HashtableTraits {
    pub const BUCKET_RECORDS_LOG2: u32 = 3;  // 8 records per bucket
    pub const HASH_RECORD_BYTES_LOG2: u32 = 3; // 8 bytes per record
    pub const HASH_BUCKET_BYTES_LOG2: u32 = 6; // 64 bytes per bucket
    pub const MAX_WRAP_BYTES_LOG2: u32 = 16; // Maximum wrap size
}

impl Bucket {
    /// Create a new empty bucket
    pub fn new() -> Self {
        Self {
            records: [HashRecord::empty(); 8],
        }
    }

    /// Get the number of records per bucket
    pub const fn records_per_bucket() -> usize {
        1 << HashtableTraits::BUCKET_RECORDS_LOG2
    }

    /// Get the size of a bucket in bytes
    pub const fn bytes_per_bucket() -> usize {
        1 << HashtableTraits::HASH_BUCKET_BYTES_LOG2
    }

    /// Get a record by index (0-7)
    pub fn get(&self, index: usize) -> Option<&HashRecord> {
        self.records.get(index)
    }

    /// Get a mutable record by index (0-7)
    pub fn get_mut(&mut self, index: usize) -> Option<&mut HashRecord> {
        self.records.get_mut(index)
    }

    /// Set a record at the given index
    pub fn set(&mut self, index: usize, record: HashRecord) -> Result<(), &'static str> {
        if index >= 8 {
            return Err("Bucket index out of range (0-7)");
        }
        self.records[index] = record;
        Ok(())
    }

    /// Find the first empty slot in the bucket
    pub fn find_empty_slot(&self) -> Option<usize> {
        self.records.iter()
            .position(|record| record.is_empty())
    }

    /// Check if the bucket is full
    pub fn is_full(&self) -> bool {
        self.records.iter()
            .all(|record| !record.is_empty())
    }

    /// Count non-empty records in the bucket
    pub fn count_occupied(&self) -> usize {
        self.records.iter()
            .filter(|record| !record.is_empty())
            .count()
    }

    /// Iterate over all records
    pub fn iter(&self) -> std::slice::Iter<HashRecord> {
        self.records.iter()
    }

    /// Iterate over all records mutably
    pub fn iter_mut(&mut self) -> std::slice::IterMut<HashRecord> {
        self.records.iter_mut()
    }

    /// Convert to raw bytes for serialization
    pub fn to_bytes(&self) -> [u8; 64] {
        let mut bytes = [0u8; 64];
        for (i, record) in self.records.iter().enumerate() {
            let record_bytes = record.raw().to_le_bytes();
            bytes[i * 8..(i + 1) * 8].copy_from_slice(&record_bytes);
        }
        bytes
    }

    /// Create from raw bytes (for deserialization)
    pub fn from_bytes(bytes: &[u8; 64]) -> Self {
        let mut records = [HashRecord::empty(); 8];
        for i in 0..8 {
            let mut record_bytes = [0u8; 8];
            record_bytes.copy_from_slice(&bytes[i * 8..(i + 1) * 8]);
            let raw_value = u64::from_le_bytes(record_bytes);
            records[i] = HashRecord::from_raw(raw_value);
        }
        Self { records }
    }
}

impl Default for Bucket {
    fn default() -> Self {
        Self::new()
    }
}

impl std::ops::Index<usize> for Bucket {
    type Output = HashRecord;

    fn index(&self, index: usize) -> &Self::Output {
        &self.records[index]
    }
}

impl std::ops::IndexMut<usize> for Bucket {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.records[index]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bucket_creation() {
        let bucket = Bucket::new();
        assert_eq!(bucket.count_occupied(), 0);
        assert!(!bucket.is_full());
        
        // All records should be empty
        for i in 0..8 {
            assert!(bucket[i].is_empty());
        }
    }

    #[test]
    fn test_bucket_operations() {
        let mut bucket = Bucket::new();
        
        // Add a HIT record
        let hit_record = HashRecord::hit(1, 0x12345, false, false, false, 1000);
        bucket.set(0, hit_record).unwrap();
        
        assert_eq!(bucket.count_occupied(), 1);
        assert_eq!(bucket.find_empty_slot(), Some(1));
        
        // Fill the bucket
        for i in 1..8 {
            let record = HashRecord::hit(i as u8, 0x10000 + i as u32, false, false, false, i as u32 * 100);
            bucket.set(i, record).unwrap();
        }
        
        assert_eq!(bucket.count_occupied(), 8);
        assert!(bucket.is_full());
        assert_eq!(bucket.find_empty_slot(), None);
    }

    #[test]
    fn test_bucket_serialization() {
        let mut bucket = Bucket::new();
        
        // Add some records
        bucket.set(0, HashRecord::hit(1, 0x12345, false, false, false, 1000)).unwrap();
        bucket.set(1, HashRecord::hifreq(2, 0x54321, true, false, false, false, 500)).unwrap();
        
        // Serialize to bytes
        let bytes = bucket.to_bytes();
        
        // Deserialize back
        let bucket2 = Bucket::from_bytes(&bytes);
        
        // Check that records match
        assert_eq!(bucket[0].raw(), bucket2[0].raw());
        assert_eq!(bucket[1].raw(), bucket2[1].raw());
        
        // Check empty records are preserved
        for i in 2..8 {
            assert!(bucket2[i].is_empty());
        }
    }

    #[test]
    fn test_bucket_constants() {
        assert_eq!(Bucket::records_per_bucket(), 8);
        assert_eq!(Bucket::bytes_per_bucket(), 64);
        assert_eq!(1 << HashtableTraits::BUCKET_RECORDS_LOG2, 8);
        assert_eq!(1 << HashtableTraits::HASH_BUCKET_BYTES_LOG2, 64);
    }
}