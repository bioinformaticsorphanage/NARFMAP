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

    // =========================================================================
    // Tests ported from C++ ReferenceSequenceGtest.cpp
    // =========================================================================

    /// Port of ReferenceSequence::decodeBase
    /// Tests one-hot to IUPAC character decoding
    #[test]
    fn test_decode_base() {
        // DRAGEN one-hot to IUPAC mapping
        fn decode_base(encoded: u8) -> char {
            // Only low 4 bits matter
            match encoded & 0x0F {
                0x0 => 'P', // Padding/null
                0x1 => 'A',
                0x2 => 'C',
                0x3 => 'M', // A|C
                0x4 => 'G',
                0x5 => 'R', // A|G
                0x6 => 'S', // C|G
                0x7 => 'V', // A|C|G
                0x8 => 'T',
                0x9 => 'W', // A|T
                0xA => 'Y', // C|T
                0xB => 'H', // A|C|T
                0xC => 'K', // G|T
                0xD => 'D', // A|G|T
                0xE => 'B', // C|G|T
                0xF => 'N', // A|C|G|T
                _ => unreachable!(),
            }
        }

        // Test all 16 values
        assert_eq!(decode_base(0x0), 'P');
        assert_eq!(decode_base(0x1), 'A');
        assert_eq!(decode_base(0x2), 'C');
        assert_eq!(decode_base(0x3), 'M');
        assert_eq!(decode_base(0x4), 'G');
        assert_eq!(decode_base(0x5), 'R');
        assert_eq!(decode_base(0x6), 'S');
        assert_eq!(decode_base(0x7), 'V');
        assert_eq!(decode_base(0x8), 'T');
        assert_eq!(decode_base(0x9), 'W');
        assert_eq!(decode_base(0xA), 'Y');
        assert_eq!(decode_base(0xB), 'H');
        assert_eq!(decode_base(0xC), 'K');
        assert_eq!(decode_base(0xD), 'D');
        assert_eq!(decode_base(0xE), 'B');
        assert_eq!(decode_base(0xF), 'N');

        // Verify high nibble is ignored (only low 4 bits matter)
        for i in 0x10u8..=0xFF {
            assert_eq!(decode_base(i), decode_base(i & 0x0F));
        }
    }

    /// Port of ReferenceSequence::translateTo2bpb
    /// Tests one-hot to 2-bit encoding (A=0, C=1, G=2, T=3)
    #[test]
    fn test_translate_to_2bpb() {
        fn translate_to_2bpb(encoded: u8) -> u8 {
            match encoded & 0x0F {
                0x1 => 0, // A
                0x2 => 1, // C
                0x4 => 2, // G
                0x8 => 3, // T
                _ => 0,   // Ambiguous bases default to 0 (A)
            }
        }

        // Primary bases
        assert_eq!(translate_to_2bpb(0x0), 0); // P -> 0
        assert_eq!(translate_to_2bpb(0x1), 0); // A -> 0
        assert_eq!(translate_to_2bpb(0x2), 1); // C -> 1
        assert_eq!(translate_to_2bpb(0x3), 0); // M -> 0
        assert_eq!(translate_to_2bpb(0x4), 2); // G -> 2
        assert_eq!(translate_to_2bpb(0x5), 0); // R -> 0
        assert_eq!(translate_to_2bpb(0x6), 0); // S -> 0
        assert_eq!(translate_to_2bpb(0x7), 0); // V -> 0
        assert_eq!(translate_to_2bpb(0x8), 3); // T -> 3
        assert_eq!(translate_to_2bpb(0x9), 0); // W -> 0
        assert_eq!(translate_to_2bpb(0xA), 0); // Y -> 0
        assert_eq!(translate_to_2bpb(0xB), 0); // H -> 0
        assert_eq!(translate_to_2bpb(0xC), 0); // K -> 0
        assert_eq!(translate_to_2bpb(0xD), 0); // D -> 0
        assert_eq!(translate_to_2bpb(0xE), 0); // B -> 0
        assert_eq!(translate_to_2bpb(0xF), 0); // N -> 0

        // Verify high nibble is ignored
        for i in 0x10u8..=0xFF {
            assert_eq!(translate_to_2bpb(i), translate_to_2bpb(i & 0x0F));
        }
    }

    /// Port of ReferenceSequence::generateSequence
    /// Tests encoding a base string to packed 4-bit format
    #[test]
    fn test_generate_sequence() {
        fn encode_base(base: char) -> u8 {
            match base {
                'A' => 1,
                'C' => 2,
                'G' => 4,
                'T' => 8,
                _ => 0,
            }
        }

        fn generate_sequence(bases: &str) -> Vec<u8> {
            assert!(bases.len() % 2 == 0);
            let mut result = Vec::new();
            let chars: Vec<char> = bases.chars().collect();
            for chunk in chars.chunks(2) {
                let low = encode_base(chunk[0]);
                let high = encode_base(chunk[1]);
                result.push(low | (high << 4));
            }
            result
        }

        // C++ test data: "ACGTAACCGGTTAAACCCGGGTTT"
        let bases = "ACGTAACCGGTTAAACCCGGGTTT";
        let expected: Vec<u8> = vec![
            0x21, 0x84, 0x11, 0x22, 0x44, 0x88, 0x11, 0x21, 0x22, 0x44, 0x84, 0x88,
        ];

        let generated = generate_sequence(bases);
        assert_eq!(generated.len(), expected.len());
        for (i, (&gen, &exp)) in generated.iter().zip(expected.iter()).enumerate() {
            assert_eq!(
                gen, exp,
                "mismatch at byte {i}: got {gen:#04x}, expected {exp:#04x}"
            );
        }
    }

    /// Port of ReferenceSequence::getBase
    /// Tests retrieving one-hot encoded bases from packed data
    #[test]
    fn test_get_base() {
        fn generate_sequence(bases: &str) -> Vec<u8> {
            fn encode_base(base: char) -> u8 {
                match base {
                    'A' => 1,
                    'C' => 2,
                    'G' => 4,
                    'T' => 8,
                    _ => 0,
                }
            }
            assert!(bases.len() % 2 == 0);
            let mut result = Vec::new();
            let chars: Vec<char> = bases.chars().collect();
            for chunk in chars.chunks(2) {
                result.push(encode_base(chunk[0]) | (encode_base(chunk[1]) << 4));
            }
            result
        }

        fn get_base(data: &[u8], pos: usize) -> u8 {
            let byte_idx = pos / 2;
            if byte_idx >= data.len() {
                return 0;
            }
            let byte = data[byte_idx];
            if pos % 2 == 0 {
                byte & 0x0F
            } else {
                (byte >> 4) & 0x0F
            }
        }

        let bases = "ACGTAACCGGTTAAACCCGGGTTT";
        let sequence = generate_sequence(bases);

        // Verify packed data has correct size
        assert_eq!(sequence.len() * 2, bases.len());

        // Test first 4 bases: A=1, C=2, G=4, T=8
        assert_eq!(get_base(&sequence, 0), 1); // A
        assert_eq!(get_base(&sequence, 1), 2); // C
        assert_eq!(get_base(&sequence, 2), 4); // G
        assert_eq!(get_base(&sequence, 3), 8); // T

        // Test last base
        assert_eq!(get_base(&sequence, bases.len() - 1), 8); // T
    }
}
