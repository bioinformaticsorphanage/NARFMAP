//! MAPQ (Mapping Quality) calculation
//!
//! Ported from DRAGEN's Mapq.hpp

/// Maximum MAPQ value
pub const MAPQ_MAX: i32 = 60;

/// Hardware MAPQ max (for debugging)
#[allow(dead_code)]
pub const HW_MAPQ_MAX: i32 = 250;

/// MAPQ coefficient (matches DRAGEN: 38912 >> 8 = 152)
const MAPQ_COEFF: i32 = 152;
const MAPQ_COEFF_I: i32 = 38912;

/// Lookup table for log2 approximation (7-bit fractional precision)
/// Matches DRAGEN hardware implementation exactly
#[allow(clippy::unreadable_literal)]
const LOG2_APPROX_TABLE: [i32; 128] = [
    0b0000000, 0b0000001, 0b0000011, 0b0000100, 0b0000110, 0b0000111, 0b0001000, 0b0001010,
    0b0001011, 0b0001101, 0b0001110, 0b0001111, 0b0010001, 0b0010010, 0b0010011, 0b0010100,
    0b0010110, 0b0010111, 0b0011000, 0b0011010, 0b0011011, 0b0011100, 0b0011101, 0b0011111,
    0b0100000, 0b0100001, 0b0100010, 0b0100011, 0b0100101, 0b0100110, 0b0100111, 0b0101000,
    0b0101001, 0b0101010, 0b0101100, 0b0101101, 0b0101110, 0b0101111, 0b0110000, 0b0110001,
    0b0110010, 0b0110011, 0b0110100, 0b0110101, 0b0110111, 0b0111000, 0b0111001, 0b0111010,
    0b0111011, 0b0111100, 0b0111101, 0b0111110, 0b0111111, 0b1000000, 0b1000001, 0b1000010,
    0b1000011, 0b1000100, 0b1000101, 0b1000110, 0b1000111, 0b1001000, 0b1001001, 0b1001010,
    0b1001011, 0b1001100, 0b1001101, 0b1001110, 0b1001111, 0b1010000, 0b1010001, 0b1010001,
    0b1010010, 0b1010011, 0b1010100, 0b1010101, 0b1010110, 0b1010111, 0b1011000, 0b1011001,
    0b1011010, 0b1011011, 0b1011011, 0b1011100, 0b1011101, 0b1011110, 0b1011111, 0b1100000,
    0b1100001, 0b1100001, 0b1100010, 0b1100011, 0b1100100, 0b1100101, 0b1100110, 0b1100111,
    0b1100111, 0b1101000, 0b1101001, 0b1101010, 0b1101011, 0b1101011, 0b1101100, 0b1101101,
    0b1101110, 0b1101111, 0b1101111, 0b1110000, 0b1110001, 0b1110010, 0b1110011, 0b1110011,
    0b1110100, 0b1110101, 0b1110110, 0b1110110, 0b1110111, 0b1111000, 0b1111001, 0b1111001,
    0b1111010, 0b1111011, 0b1111100, 0b1111100, 0b1111101, 0b1111110, 0b1111111, 0b1111111,
];

/// Log2 approximation matching DRAGEN hardware
/// Returns fixed-point result with 7 fractional bits
pub fn log2_approx(x: i32) -> i32 {
    if x <= 0 {
        return 0;
    }

    // Find integer portion of log (position of MSB)
    let log_int = 31 - x.leading_zeros() as i32;

    // Normalize to [1,2) and get 7 fractional bits
    let norm = ((x << 7) >> log_int) as usize;

    // Lookup fractional part
    let log_frac = LOG2_APPROX_TABLE[norm & 0x7f];

    // Combine integer and fractional portions
    (log_int << 7) + log_frac
}

/// MAPQ coefficient scaled by SNP cost
fn mapq_coeff_scaled_i(snp_cost: i32) -> i32 {
    MAPQ_COEFF_I * 5 / snp_cost
}

/// Alignment score to MAPQ scaling factor
fn aln2mapq(snp_cost: i32, read_len: i32) -> i32 {
    let log2_length = log2_approx(read_len);
    let coeff = (MAPQ_COEFF as f64) * (5.0 / snp_cost as f64);
    let divisor = ((log2_length * log2_length) >> 7) as f64;
    if divisor == 0.0 {
        return i32::MAX;
    }
    (coeff / divisor * (1 << 20) as f64) as i32
}

/// Compute MAPQ from alignment scores
///
/// # Arguments
/// * `snp_cost` - SNP mismatch penalty (typically 5)
/// * `as_score` - Alignment score of best alignment
/// * `xs_score` - Alignment score of second-best alignment (0 if none)
/// * `read_len` - Read length
///
/// # Returns
/// MAPQ value capped at MAPQ_MAX (60)
pub fn compute_mapq(snp_cost: i32, as_score: i32, xs_score: i32, read_len: i32) -> i32 {
    let s1 = as_score;
    let s2 = xs_score;

    let a2m_scale = aln2mapq(snp_cost, read_len);
    let mapq = ((s1 - s2) as i64 * a2m_scale as i64) >> 13;

    mapq.min(MAPQ_MAX as i64) as i32
}

/// MAPQ to alignment score conversion
#[allow(dead_code)]
pub fn mapq2aln(snp_cost: i32, read_length: i32) -> i32 {
    let log2_length = log2_approx(read_length);
    log2_length * log2_length / (mapq_coeff_scaled_i(snp_cost) >> 4)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_log2_approx_powers_of_two() {
        // log2(1) = 0.0 -> 0 << 7 = 0
        assert_eq!(log2_approx(1), 0);
        // log2(2) = 1.0 -> 1 << 7 = 128
        assert_eq!(log2_approx(2), 128);
        // log2(4) = 2.0 -> 2 << 7 = 256
        assert_eq!(log2_approx(4), 256);
        // log2(8) = 3.0 -> 3 << 7 = 384
        assert_eq!(log2_approx(8), 384);
    }

    #[test]
    fn test_log2_approx_typical_read_lengths() {
        // log2(100) ≈ 6.64 -> 6 << 7 + frac ≈ 768 + 82 = 850
        let log100 = log2_approx(100);
        assert!(log100 > 768 && log100 < 896, "log2(100) = {log100}");

        // log2(150) ≈ 7.23 -> 7 << 7 + frac ≈ 896 + 30 = 926
        let log150 = log2_approx(150);
        assert!(log150 > 896 && log150 < 1024, "log2(150) = {log150}");
    }

    #[test]
    fn test_mapq_unique_alignment() {
        // Unique alignment (no second best) should get high MAPQ
        // as=100, xs=0, read_len=100, snp_cost=5
        let mapq = compute_mapq(5, 100, 0, 100);
        assert_eq!(mapq, MAPQ_MAX, "unique alignment should get MAPQ_MAX");
    }

    #[test]
    fn test_mapq_identical_scores() {
        // Identical scores should give MAPQ = 0
        let mapq = compute_mapq(5, 100, 100, 100);
        assert_eq!(mapq, 0, "identical scores should give MAPQ 0");
    }

    #[test]
    fn test_mapq_small_difference() {
        // Small score difference should give low but non-zero MAPQ
        let mapq = compute_mapq(5, 100, 95, 100);
        assert!(mapq > 0 && mapq < MAPQ_MAX, "small diff MAPQ = {mapq}");
    }

    #[test]
    fn test_mapq_zero_second_score() {
        // No second alignment (xs=0) should give max MAPQ
        let mapq = compute_mapq(5, 50, 0, 100);
        assert_eq!(mapq, MAPQ_MAX);
    }
}
