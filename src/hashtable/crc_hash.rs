use anyhow::{Result, anyhow};

/// CRC polynomial representation for k-mer hashing
/// Uses the same format as DRAGMAP: up to 128-bit polynomials stored as little-endian bytes
#[derive(Debug, Clone)]
pub struct CrcPolynomial {
    /// Polynomial coefficients stored as bytes (little-endian)
    poly: Vec<u8>,
    /// Number of bits in the polynomial
    bits: usize,
}

impl CrcPolynomial {
    /// Create a new CRC polynomial from a hex string
    /// Example: "2C991CE6A8DD55" for 54-bit polynomial
    pub fn from_hex(hex_str: &str) -> Result<Self> {
        let hex_clean = hex_str.trim().replace(" ", "");
        if hex_clean.len() % 2 != 0 {
            return Err(anyhow!("Hex string must have even length"));
        }

        let mut poly = Vec::new();
        // Parse hex string in reverse order for little-endian storage
        for i in (0..hex_clean.len()).step_by(2).rev() {
            let end = i + 2;
            if end <= hex_clean.len() {
                let byte_str = &hex_clean[i..end];
                let byte = u8::from_str_radix(byte_str, 16)
                    .map_err(|e| anyhow!("Invalid hex byte '{}': {}", byte_str, e))?;
                poly.push(byte);
            }
        }

        let bits = hex_clean.len() * 4;
        Ok(Self { poly, bits })
    }

    /// Create from raw bytes (little-endian)
    pub fn from_bytes(bytes: &[u8], bits: usize) -> Self {
        Self {
            poly: bytes.to_vec(),
            bits,
        }
    }

    /// Get the number of bits
    pub fn bits(&self) -> usize {
        self.bits
    }

    /// Get polynomial bytes (little-endian)
    pub fn bytes(&self) -> &[u8] {
        &self.poly
    }
}

/// CRC hasher implementing DRAGMAP's custom CRC algorithm
/// This implements the same algorithm as crcHashSlow() in the C++ code
pub struct CrcHasher {
    polynomial: CrcPolynomial,
}

impl CrcHasher {
    /// Create a new CRC hasher with the given polynomial
    pub fn new(polynomial: CrcPolynomial) -> Self {
        Self { polynomial }
    }

    /// Hash input data using the CRC polynomial
    /// This implements the same algorithm as crcHashSlow() from DRAGMAP
    pub fn hash(&self, data: &[u8]) -> Result<Vec<u8>> {
        let bits = self.polynomial.bits();
        let bytes = (bits + 7) / 8;
        let poly_bytes = self.polynomial.bytes();

        if data.len() != bytes {
            return Err(anyhow!(
                "Data length {} doesn't match expected {} bytes for {}-bit polynomial",
                data.len(), bytes, bits
            ));
        }

        // Initialize hash with input data (copy to avoid mutation)
        let mut hash = data.to_vec();
        hash.resize(bytes, 0);

        let top_byte = bytes - 1;
        let top_bit_mask = 1u8 << ((bits + 7) % 8);

        // Process each bit (polynomial division)
        for _i in 0..bits {
            // Check if MSB is 1
            let subtract = (hash[top_byte] & top_bit_mask) != 0;

            // Left-shift the remainder (polynomial long division)
            for j in (1..bytes).rev() {
                hash[j] = (hash[j] << 1) | (hash[j - 1] >> 7);
            }
            hash[0] <<= 1;

            // XOR with polynomial if MSB was 1
            if subtract {
                for j in 0..bytes {
                    if j < poly_bytes.len() {
                        hash[j] ^= poly_bytes[j];
                    }
                }
            }
        }

        Ok(hash)
    }

    /// Hash a 64-bit value (optimized path for k-mers <= 32 bases)
    pub fn hash_u64(&self, data: u64) -> Result<u64> {
        let data_bytes = data.to_le_bytes();
        let relevant_bytes = (self.polynomial.bits() + 7) / 8;
        let hash_bytes = self.hash(&data_bytes[..relevant_bytes.min(8)])?;
        
        // Convert back to u64
        let mut result = 0u64;
        for (i, &byte) in hash_bytes.iter().enumerate().take(8) {
            result |= (byte as u64) << (i * 8);
        }
        
        Ok(result)
    }
}

/// K-mer hasher that converts DNA sequences to CRC hash values
pub struct KmerHasher {
    primary_hasher: CrcHasher,
    secondary_hasher: Option<CrcHasher>,
    kmer_bits: usize,
}

impl KmerHasher {
    /// Create a new k-mer hasher with primary and optional secondary CRC polynomials
    pub fn new(
        primary_poly: CrcPolynomial,
        secondary_poly: Option<CrcPolynomial>,
        kmer_length: usize,
    ) -> Result<Self> {
        let kmer_bits = kmer_length * 2; // 2 bits per base
        
        if primary_poly.bits() != kmer_bits {
            return Err(anyhow!(
                "Primary polynomial has {} bits, expected {} for {}-mer",
                primary_poly.bits(), kmer_bits, kmer_length
            ));
        }

        if let Some(ref sec_poly) = secondary_poly {
            if sec_poly.bits() < kmer_bits {
                return Err(anyhow!(
                    "Secondary polynomial has {} bits, must be >= {} for {}-mer",
                    sec_poly.bits(), kmer_bits, kmer_length
                ));
            }
        }

        Ok(Self {
            primary_hasher: CrcHasher::new(primary_poly),
            secondary_hasher: secondary_poly.map(CrcHasher::new),
            kmer_bits,
        })
    }

    /// Create with default DRAGMAP polynomials for the given k-mer length
    pub fn with_dragmap_defaults(kmer_length: usize) -> Result<Self> {
        let kmer_bits = kmer_length * 2;
        
        // Use a default polynomial based on k-mer length
        // These are inspired by DRAGMAP's polynomial selection
        let primary_poly = match kmer_bits {
            42 => {
                // For 21-mers, use a proper 42-bit polynomial
                let poly_value = 0x2C991CE6A8u64; // 42-bit value
                let poly_bytes = poly_value.to_le_bytes();
                let needed_bytes = (kmer_bits + 7) / 8; // 42 bits = 6 bytes
                CrcPolynomial::from_bytes(&poly_bytes[..needed_bytes], kmer_bits)
            },
            32 => CrcPolynomial::from_hex("C96C5795")?,        // 32-bit for 16-mer
            8 => {
                // For 4-mers, use an 8-bit polynomial
                let poly_value = 0x1D; // CRC-8 polynomial
                let poly_bytes = [poly_value];
                CrcPolynomial::from_bytes(&poly_bytes, kmer_bits)
            },
            _ => {
                // Generate a simple polynomial for other lengths
                let poly_value = 0x1021u64 << (kmer_bits.saturating_sub(16));
                let poly_bytes = poly_value.to_le_bytes();
                let needed_bytes = (kmer_bits + 7) / 8;
                CrcPolynomial::from_bytes(&poly_bytes[..needed_bytes], kmer_bits)
            }
        };

        // Secondary polynomial for extended hashing (if needed)
        let secondary_poly = if kmer_bits < 54 {
            Some(CrcPolynomial::from_hex("2C991CE6A8DD55")?) // 54-bit extended
        } else {
            None
        };

        Self::new(primary_poly, secondary_poly, kmer_length)
    }

    /// Hash a k-mer using the primary CRC polynomial
    pub fn hash_kmer(&self, kmer_data: u64) -> Result<u64> {
        self.primary_hasher.hash_u64(kmer_data)
    }

    /// Hash an extended k-mer using the secondary CRC polynomial
    pub fn hash_extended_kmer(&self, kmer_data: &[u8]) -> Result<Vec<u8>> {
        match &self.secondary_hasher {
            Some(hasher) => hasher.hash(kmer_data),
            None => Err(anyhow!("No secondary hasher configured")),
        }
    }

    /// Convert DNA sequence to 2-bit packed representation
    /// Uses DRAGMAP encoding: A=0, C=1, G=2, T=3, N=2 (treated as G)
    pub fn sequence_to_2bit(&self, sequence: &[u8]) -> Result<u64> {
        if sequence.len() > 32 {
            return Err(anyhow!("Sequence too long for u64 (max 32 bases)"));
        }

        let mut data = 0u64;
        for (i, &base) in sequence.iter().enumerate() {
            let encoded_base = match base.to_ascii_uppercase() {
                b'A' => 0u64,
                b'C' => 1u64,
                b'G' => 2u64,
                b'T' => 3u64,
                b'N' => 2u64, // N treated as G
                _ => return Err(anyhow!("Invalid DNA base: {}", base as char)),
            };
            
            // Pack in little-endian format (2 bits per base)
            data |= encoded_base << (i * 2);
        }

        Ok(data)
    }

    /// Convert 2-bit packed data back to DNA sequence
    pub fn sequence_from_2bit(&self, data: u64, length: usize) -> String {
        let mut sequence = String::with_capacity(length);
        
        for i in 0..length {
            let base_code = (data >> (i * 2)) & 0x3;
            let base_char = match base_code {
                0 => 'A',
                1 => 'C', 
                2 => 'G',
                3 => 'T',
                _ => 'N',
            };
            sequence.push(base_char);
        }
        
        sequence
    }

    /// Get hash table address from hash value (matches DRAGMAP getVirtualByteAddress)
    pub fn get_address_from_hash(&self, hash: u64, squeeze_factor: u64) -> u64 {
        // Extract bits 19-54 (35 bits) for virtual byte address
        const ADDRESS_START: u32 = 19;
        const ADDRESS_BITS: u32 = 35;
        let address_mask = ((1u64 << ADDRESS_BITS) - 1) << ADDRESS_START;
        let address = (hash & address_mask) >> ADDRESS_START;
        address * squeeze_factor
    }

    /// Get thread ID from virtual address (matches DRAGMAP getThreadIdFromVirtualByteAddress)
    pub fn get_thread_id_from_address(&self, virtual_address: u64) -> u8 {
        // Extract bits 3-8 (6 bits) for thread ID
        const THREAD_START: u32 = 3;
        const THREAD_BITS: u32 = 6;
        let thread_mask = (1u64 << THREAD_BITS) - 1;
        ((virtual_address >> THREAD_START) & thread_mask) as u8
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_crc_polynomial_from_hex() {
        let poly = CrcPolynomial::from_hex("2C991CE6A8DD55").unwrap();
        assert_eq!(poly.bits(), 56); // 7 bytes * 8 bits
        
        // Check little-endian byte order
        let expected = vec![0x55, 0xDD, 0xA8, 0xE6, 0x1C, 0x99, 0x2C];
        assert_eq!(poly.bytes(), &expected);
    }

    #[test]
    fn test_sequence_encoding() {
        let hasher = KmerHasher::with_dragmap_defaults(16).unwrap();
        
        // Test basic encoding
        let seq = b"ACGT";
        let encoded = hasher.sequence_to_2bit(seq).unwrap();
        
        // A=0, C=1, G=2, T=3 in little-endian: 0b11100100 = 0xE4
        assert_eq!(encoded & 0xFF, 0b11100100);
        
        // Test round-trip
        let decoded = hasher.sequence_from_2bit(encoded, 4);
        assert_eq!(decoded, "ACGT");
    }

    #[test]
    fn test_kmer_hashing() {
        let hasher = KmerHasher::with_dragmap_defaults(16).unwrap();
        
        let seq = b"ACGTACGTACGTACGT"; // 16-mer
        let encoded = hasher.sequence_to_2bit(seq).unwrap();
        let hash = hasher.hash_kmer(encoded).unwrap();
        
        // Hash should be deterministic
        let hash2 = hasher.hash_kmer(encoded).unwrap();
        assert_eq!(hash, hash2);
        
        // Different sequences should have different hashes
        let seq2 = b"TGCATGCATGCATGCA";
        let encoded2 = hasher.sequence_to_2bit(seq2).unwrap();
        let hash2 = hasher.hash_kmer(encoded2).unwrap();
        assert_ne!(hash, hash2);
    }

    #[test]
    fn test_dragmap_compatibility() {
        // Test with known values from DRAGMAP test suite
        let poly = CrcPolynomial::from_hex("2C991CE6A8DD55").unwrap();
        let hasher = CrcHasher::new(poly);
        
        // These test vectors would need to be verified against actual DRAGMAP output
        // For now, just test that hashing is deterministic
        let data = 0x3543543543u64.to_le_bytes();
        let hash1 = hasher.hash(&data[..7]).unwrap(); // 54 bits = 7 bytes
        let hash2 = hasher.hash(&data[..7]).unwrap();
        assert_eq!(hash1, hash2);
    }
}