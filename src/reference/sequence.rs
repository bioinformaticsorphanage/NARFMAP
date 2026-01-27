//! Reference sequence loading and access

use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;

use anyhow::{Context, Result};

/// Reference sequence information
#[derive(Debug, Clone)]
pub struct RefSeqInfo {
    pub name: String,
    pub length: u64,
    pub offset: u64, // offset in binary reference file
}

/// Reference sequence data
pub struct ReferenceSequence {
    pub sequences: Vec<RefSeqInfo>,
    data: Vec<u8>,
}

impl ReferenceSequence {
    /// Load reference from reference directory
    pub fn load(ref_dir: &Path) -> Result<Self> {
        let ref_path = ref_dir.join("reference.bin");
        let mut file = File::open(&ref_path)
            .with_context(|| format!("Failed to open reference: {}", ref_path.display()))?;

        let mut data = Vec::new();
        file.read_to_end(&mut data)?;

        // For now, assume single sequence
        // TODO: Parse ref_index.bin for proper sequence info
        let sequences = vec![RefSeqInfo {
            name: "ref".to_string(),
            length: data.len() as u64 * 4, // 2 bits per base, 4 bases per byte
            offset: 0,
        }];

        Ok(Self { sequences, data })
    }

    /// Get sequence at position (returns bases as ASCII)
    pub fn get_sequence(&self, start: u64, length: usize) -> Vec<u8> {
        let mut result = Vec::with_capacity(length);

        for i in 0..length {
            let pos = start + i as u64;
            let byte_idx = (pos / 4) as usize;
            let bit_offset = ((pos % 4) * 2) as u32;

            if byte_idx >= self.data.len() {
                result.push(b'N');
                continue;
            }

            let bits = (self.data[byte_idx] >> (6 - bit_offset)) & 0x03;
            let base = match bits {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => b'N',
            };
            result.push(base);
        }

        result
    }

    /// Get total reference length
    pub fn total_length(&self) -> u64 {
        self.sequences.iter().map(|s| s.length).sum()
    }

    /// Find sequence containing position
    pub fn find_sequence(&self, pos: u64) -> Option<(&RefSeqInfo, u64)> {
        let mut offset = 0u64;
        for seq in &self.sequences {
            if pos < offset + seq.length {
                return Some((seq, pos - offset));
            }
            offset += seq.length;
        }
        None
    }
}
