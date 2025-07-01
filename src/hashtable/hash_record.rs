use std::fmt;

/// DRAGMAP-compatible hash table record (64-bit)
/// 
/// Record types and bit layouts:
/// - HIT: ThreadId[63:58] HashBits[57:35] EX[34] LF[33] RC[32] RefPos[31:0]
/// - EMPTY: 0[63:32] 0xF[31:28] 0x0[27:24] 0[23:0]
/// - HIFREQ: ThreadId[63:58] HashBits[57:35] EX[34] LF[33] RC[32] 0xF[31:28] 0x1[27:24] 0[23] AL[22] Frequency[21:0]
/// - EXTEND: ThreadId[63:58] HashBits[57:35] EX[34] LF[33] RC[32] 0xF[31:28] 0x2[27:24] RF[23] AL[22] ExtensionLength[21:18] ExtensionId[17:0]
/// - CHAIN_BEG_MASK: FilterMask[63:32] 0xF[31:28] 0x4[27:24] 0[23:18] ChainPointer[17:0]
/// - CHAIN_BEG_LIST: FilterList4[63:56] FilterList3[55:48] FilterList2[47:40] FilterList1[39:32] 0xF[31:28] 0x5[27:24] 0[23:18] ChainPointer[17:0]
/// - CHAIN_CON_MASK: FilterMask[63:32] 0xF[31:28] 0x6[27:24] 0[23:18] ChainPointer[17:0]
/// - CHAIN_CON_LIST: FilterList4[63:56] FilterList3[55:48] FilterList2[47:40] FilterList1[39:32] 0xF[31:28] 0x7[27:24] 0[23:18] ChainPointer[17:0]
/// - INTERVAL_SL: ThreadId[63:58] HashBits[57:35] EX[34] LF[33] RC[32] 0xF[31:28] 0x8[27:24] Length[23:15] Start[14:0]
/// - INTERVAL_SLE: ThreadId[63:58] HashBits[57:35] EX[34] LF[33] MSB[32] 0xF[31:28] 0x9[27:24] Exlifts[23:16] Length[15:8] Start[7:0]
/// - INTERVAL_S: ThreadId[63:58] HashBits[57:35] EX[34] LF[33] MSB[32] 0xF[31:28] 0xA[27:24] Start[23:0]
/// - INTERVAL_L: ThreadId[63:58] HashBits[57:35] EX[34] LF[33] MSB[32] 0xF[31:28] 0xB[27:24] Length[23:0]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HashRecord(u64);

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RecordType {
    Empty = 0,
    HiFreq = 1,
    Extend = 2,
    Repair = 3, // Obsolete
    ChainBegMask = 4,
    ChainBegList = 5,
    ChainConMask = 6,
    ChainConList = 7,
    IntervalSL = 8,
    IntervalSLE = 9,
    IntervalS = 0xA,
    IntervalL = 0xB,
    Hit = 16, // Default type when bits [31:28] are not all set
}

impl HashRecord {
    // Bit field constants
    const THREAD_ID_START: u32 = 58;
    const THREAD_ID_BITS: u32 = 6;
    const HASH_BITS_START: u32 = 35;
    const HASH_BITS_BITS: u32 = 23;
    const EX_FLAG: u32 = 34; // Extended seed flag
    const LF_FLAG: u32 = 33; // Last in thread flag
    const RC_FLAG: u32 = 32; // Reverse complement flag
    const RS_FLAG: u32 = 32; // Random sample flag (HIFREQ/EXTEND only)
    const REFERENCE_POSITION_START: u32 = 0;
    const REFERENCE_POSITION_BITS: u32 = 32;
    const NOT_HIT_START: u32 = 28;
    const NOT_HIT_BITS: u32 = 4;
    const OP_CODE_START: u32 = 24;
    const OP_CODE_BITS: u32 = 4;
    const AL_FLAG: u32 = 22; // Alt-aware flag
    const FREQUENCY_START: u32 = 0;
    const FREQUENCY_BITS: u32 = 22;
    const CHAIN_POINTER_START: u32 = 0;
    const CHAIN_POINTER_BITS: u32 = 18;

    /// Create an empty hash record
    pub fn empty() -> Self {
        Self(0xF0000000)
    }

    /// Create a HIT record
    pub fn hit(thread_id: u8, hash_bits: u32, is_extended: bool, is_last: bool, is_rc: bool, position: u32) -> Self {
        let mut record = 0u64;
        record |= (thread_id as u64 & 0x3F) << Self::THREAD_ID_START;
        record |= (hash_bits as u64 & 0x7FFFFF) << Self::HASH_BITS_START;
        if is_extended {
            record |= 1u64 << Self::EX_FLAG;
        }
        if is_last {
            record |= 1u64 << Self::LF_FLAG;
        }
        if is_rc {
            record |= 1u64 << Self::RC_FLAG;
        }
        record |= position as u64;
        Self(record)
    }

    /// Create a HIFREQ record
    pub fn hifreq(thread_id: u8, hash_bits: u32, is_extended: bool, is_last: bool, has_random_sample: bool, is_alt: bool, frequency: u32) -> Self {
        let mut record = 0u64;
        record |= (thread_id as u64 & 0x3F) << Self::THREAD_ID_START;
        record |= (hash_bits as u64 & 0x7FFFFF) << Self::HASH_BITS_START;
        if is_extended {
            record |= 1u64 << Self::EX_FLAG;
        }
        if is_last {
            record |= 1u64 << Self::LF_FLAG;
        }
        if has_random_sample {
            record |= 1u64 << Self::RS_FLAG;
        }
        record |= 0xF1000000u64; // Set bits for HIFREQ type
        if is_alt {
            record |= 1u64 << Self::AL_FLAG;
        }
        record |= (frequency as u64) & 0x3FFFFF;
        Self(record)
    }

    /// Get the record type
    pub fn record_type(&self) -> RecordType {
        // Check if bits [31:28] are all set
        if (self.0 >> Self::NOT_HIT_START) & 0xF == 0xF {
            // Extract op code from bits [27:24]
            let op_code = ((self.0 >> Self::OP_CODE_START) & 0xF) as u8;
            match op_code {
                0 => RecordType::Empty,
                1 => RecordType::HiFreq,
                2 => RecordType::Extend,
                3 => RecordType::Repair,
                4 => RecordType::ChainBegMask,
                5 => RecordType::ChainBegList,
                6 => RecordType::ChainConMask,
                7 => RecordType::ChainConList,
                8 => RecordType::IntervalSL,
                9 => RecordType::IntervalSLE,
                0xA => RecordType::IntervalS,
                0xB => RecordType::IntervalL,
                _ => RecordType::Hit, // Shouldn't happen
            }
        } else {
            RecordType::Hit
        }
    }

    /// Check if this is an empty record
    pub fn is_empty(&self) -> bool {
        matches!(self.record_type(), RecordType::Empty)
    }

    /// Get thread ID (bits 63:58)
    pub fn thread_id(&self) -> u8 {
        ((self.0 >> Self::THREAD_ID_START) & 0x3F) as u8
    }

    /// Get hash bits (bits 57:35)
    pub fn hash_bits(&self) -> u32 {
        ((self.0 >> Self::HASH_BITS_START) & 0x7FFFFF) as u32
    }

    /// Check if extended seed (bit 34)
    pub fn is_extended(&self) -> bool {
        (self.0 >> Self::EX_FLAG) & 1 == 1
    }

    /// Check if last in thread (bit 33)
    pub fn is_last_in_thread(&self) -> bool {
        (self.0 >> Self::LF_FLAG) & 1 == 1
    }

    /// Check if reverse complement (bit 32) - only for HIT records
    pub fn is_reverse_complement(&self) -> bool {
        matches!(self.record_type(), RecordType::Hit) && ((self.0 >> Self::RC_FLAG) & 1 == 1)
    }

    /// Get reference position for HIT records
    pub fn reference_position(&self) -> Option<u32> {
        match self.record_type() {
            RecordType::Hit => Some((self.0 & 0xFFFFFFFF) as u32),
            _ => None,
        }
    }

    /// Get frequency for HIFREQ records
    pub fn frequency(&self) -> Option<u32> {
        match self.record_type() {
            RecordType::HiFreq => Some((self.0 & 0x3FFFFF) as u32),
            _ => None,
        }
    }

    /// Get chain pointer for CHAIN records
    pub fn chain_pointer(&self) -> Option<u32> {
        match self.record_type() {
            RecordType::ChainBegMask | RecordType::ChainBegList | 
            RecordType::ChainConMask | RecordType::ChainConList => {
                Some((self.0 & 0x3FFFF) as u32)
            }
            _ => None,
        }
    }

    /// Get the raw 64-bit value
    pub fn raw(&self) -> u64 {
        self.0
    }

    /// Create from raw 64-bit value
    pub fn from_raw(value: u64) -> Self {
        Self(value)
    }
}

impl fmt::Display for HashRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.record_type() {
            RecordType::Empty => write!(f, "EMPTY"),
            RecordType::Hit => {
                write!(f, "HIT(tid:{}, hash:{:06x}, ext:{}, last:{}, rc:{}, pos:{})",
                    self.thread_id(),
                    self.hash_bits(),
                    self.is_extended(),
                    self.is_last_in_thread(),
                    self.is_reverse_complement(),
                    self.reference_position().unwrap_or(0)
                )
            }
            RecordType::HiFreq => {
                write!(f, "HIFREQ(tid:{}, hash:{:06x}, freq:{})",
                    self.thread_id(),
                    self.hash_bits(),
                    self.frequency().unwrap_or(0)
                )
            }
            _ => write!(f, "{:?}({:016x})", self.record_type(), self.0),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_record() {
        let record = HashRecord::empty();
        assert_eq!(record.record_type(), RecordType::Empty);
        assert!(record.is_empty());
        assert_eq!(record.raw(), 0xF0000000);
    }

    #[test]
    fn test_hit_record() {
        let record = HashRecord::hit(
            0x15,    // thread_id
            0x12345, // hash_bits
            true,    // is_extended
            false,   // is_last
            true,    // is_rc
            0x1000   // position
        );
        
        assert_eq!(record.record_type(), RecordType::Hit);
        assert_eq!(record.thread_id(), 0x15);
        assert_eq!(record.hash_bits(), 0x12345);
        assert!(record.is_extended());
        assert!(!record.is_last_in_thread());
        assert!(record.is_reverse_complement());
        assert_eq!(record.reference_position(), Some(0x1000));
    }

    #[test]
    fn test_hifreq_record() {
        let record = HashRecord::hifreq(
            0x3F,     // thread_id
            0x7FFFFF, // hash_bits (max value)
            false,    // is_extended
            true,     // is_last
            false,    // has_random_sample
            true,     // is_alt
            1000      // frequency
        );
        
        assert_eq!(record.record_type(), RecordType::HiFreq);
        assert_eq!(record.thread_id(), 0x3F);
        assert_eq!(record.hash_bits(), 0x7FFFFF);
        assert!(!record.is_extended());
        assert!(record.is_last_in_thread());
        assert_eq!(record.frequency(), Some(1000));
    }

    #[test]
    fn test_record_type_detection() {
        // Test that 0x00000000F0000000 is detected as EMPTY
        let empty = HashRecord::from_raw(0x00000000F0000000);
        assert_eq!(empty.record_type(), RecordType::Empty);

        // Test that 0x00000000F1000000 is detected as HIFREQ
        let hifreq = HashRecord::from_raw(0x00000000F1000000);
        assert_eq!(hifreq.record_type(), RecordType::HiFreq);

        // Test that a normal value is detected as HIT
        let hit = HashRecord::from_raw(0x123456789ABCDEF0);
        assert_eq!(hit.record_type(), RecordType::Hit);
    }
}