//! Reference sequence loading and access

use std::fs::File;
use std::io::Read;
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
        // DRAGEN uses 4 bits per base (2 bases per byte)
        let sequences = vec![RefSeqInfo {
            name: "ref".to_string(),
            length: data.len() as u64 * 2, // 4 bits per base, 2 bases per byte
            offset: 0,
        }];

        Ok(Self { sequences, data })
    }

    /// Get sequence at position (returns bases as ASCII)
    /// DRAGEN format: 4 bits per base, 2 bases per byte
    /// Low nibble = even position, high nibble = odd position
    /// Base encoding: A=1, C=2, G=4, T=8 (one-hot), 0=N
    pub fn get_sequence(&self, start: u64, length: usize) -> Vec<u8> {
        let mut result = Vec::with_capacity(length);

        for i in 0..length {
            let pos = start + i as u64;
            let byte_idx = (pos / 2) as usize;
            let is_odd = (pos % 2) == 1;

            if byte_idx >= self.data.len() {
                result.push(b'N');
                continue;
            }

            let byte = self.data[byte_idx];
            // Low nibble = even positions, high nibble = odd positions
            let nibble = if is_odd {
                (byte >> 4) & 0x0F
            } else {
                byte & 0x0F
            };

            // DRAGEN one-hot encoding: A=1, C=2, G=4, T=8
            let base = match nibble {
                1 => b'A',
                2 => b'C',
                4 => b'G',
                8 => b'T',
                _ => b'N', // 0 or ambiguous
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

#[cfg(test)]
mod tests {
    use super::*;

    fn reference_with_data(data: Vec<u8>) -> ReferenceSequence {
        ReferenceSequence {
            sequences: vec![RefSeqInfo {
                name: "ref".to_string(),
                length: data.len() as u64 * 2,
                offset: 0,
            }],
            data,
        }
    }

    #[test]
    fn decodes_even_and_odd_nibbles() {
        let reference = reference_with_data(vec![0x21, 0x84]);
        let seq = reference.get_sequence(0, 4);

        assert_eq!(seq, b"ACGT");
    }

    #[test]
    fn returns_n_for_out_of_range_positions() {
        let reference = reference_with_data(vec![0x21]);
        let seq = reference.get_sequence(0, 4);

        assert_eq!(seq, b"ACNN");
    }

    #[test]
    fn find_sequence_returns_relative_position() {
        let reference = ReferenceSequence {
            sequences: vec![
                RefSeqInfo {
                    name: "chr1".to_string(),
                    length: 5,
                    offset: 0,
                },
                RefSeqInfo {
                    name: "chr2".to_string(),
                    length: 3,
                    offset: 5,
                },
            ],
            data: vec![0; 4],
        };

        let (seq, pos) = reference.find_sequence(0).unwrap();
        assert_eq!(seq.name, "chr1");
        assert_eq!(pos, 0);

        let (seq, pos) = reference.find_sequence(5).unwrap();
        assert_eq!(seq.name, "chr2");
        assert_eq!(pos, 0);

        let (seq, pos) = reference.find_sequence(7).unwrap();
        assert_eq!(seq.name, "chr2");
        assert_eq!(pos, 2);

        assert!(reference.find_sequence(8).is_none());
    }
}
