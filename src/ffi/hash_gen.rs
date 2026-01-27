//! FFI bindings to hash table generation C library
//!
//! Based on gen_hash_table.h from thirdparty/dragen

#![allow(unsafe_code)] // FFI module requires unsafe

use std::os::raw::{c_char, c_double, c_int};

/// Hash table type for generation
#[repr(C)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)] // FFI enum - variants used by C code
pub enum HashTableType {
    Normal = 0,
    MethylGToA = 1,
    MethylCToT = 2,
    MethylCombined = 3,
    Anchored = 4,
}

/// Hash table header structure (512 bytes, packed)
/// This is written to hash_table.cfg.bin
#[repr(C, packed)]
#[derive(Debug, Clone, Copy)]
pub struct HashTableHeader {
    pub hash_table_version: u32,
    pub hash_table_bytes: u64,
    pub pri_seed_bases: u32,
    pub max_seed_bases: u32,
    pub max_ext_increment: u32,
    pub ref_seed_interval: c_double,
    pub table_addr_bits: u32,
    pub table_size_64ths: u32,
    pub max_seed_freq: u32,
    pub pri_max_seed_freq: u32,
    pub max_seed_freq_len: u32,
    pub target_seed_freq: c_double,
    pub thinning_freq_cap: c_double,
    pub thinning_period: u32,
    pub pri_crc_bits: u32,
    pub sec_crc_bits: u32,
    pub seed_len_cost: c_double,
    pub seed_freq_cost: c_double,
    pub extension_cost: c_double,
    pub ext_step_cost: c_double,
    pub repair_strategy: u32,
    pub min_repair_prob: c_double,
    pub anchor_bin_bits: u32,
    pub hi_freq_rand_hit: u32,
    pub ext_rand_hit_freq: u32,
    pub pri_crc_poly: [u8; 8],
    pub sec_crc_poly: [u8; 8],
    pub ref_seq_len: u64,
    pub ref_len_raw: u64,
    pub ref_len_not_n: u64,
    pub digest: u32,
    pub num_ref_seqs: u32,
    pub digest_type: u32,
    pub ref_digest: u32,
    pub ref_index_digest: u32,
    pub hash_digest: u32,
    pub liftover_digest: u32,
    pub ref_alt_seed: u32,
    pub ref_alt_start: u64,
    pub ext_tab_recs: u32,
    pub ext_tab_digest: u32,
    pub ext_rec_cost: c_double,
    pub min_freq_to_extend: u32,
    pub max_mult_base_seeds: u32,
    pub pop_snps_digest: u32,
    pub lift_match_seed_int: u32,
    pub padding: [u8; 264],
}

/// Hash table sequence info
#[repr(C, packed)]
#[derive(Debug, Clone, Copy)]
pub struct HashTableSeq {
    pub seq_start: u64,
    pub beg_trim: u32,
    pub end_trim: u32,
    pub seq_len: u32,
}

/// Hash table configuration structure
#[repr(C)]
pub struct HashTableConfig {
    pub hdr: *mut HashTableHeader,
    pub ref_seq: *mut *mut HashTableSeq,
    pub seq_name: *mut *mut c_char,

    pub max_threads: c_int,
    pub max_gb: c_int,
    pub write_hash_file: c_int,
    pub write_comp_file: c_int,
    pub size_str: *const c_char,
    pub mem_size_str: *const c_char,
    pub sj_size_str: *const c_char,
    pub methylated_conv: u32,
    pub ext_table_alloc: u32,
    pub pri_poly_index: c_int,
    pub sec_poly_index: c_int,

    pub ref_input: *mut c_char,
    pub alt_liftover: *mut c_char,
    pub config_fname: *mut c_char,
    pub config_bin_fname: *mut c_char,
    pub hash_fname: *mut c_char,
    pub comp_fname: *mut c_char,
    pub ext_tab_fname: *mut c_char,
    pub ref_output: *mut c_char,
    pub ref_idx_fname: *mut c_char,
    pub rep_mask_fname: *mut c_char,
    pub str_fname: *mut c_char,
    pub stats_fname: *mut c_char,
    pub decoy_fname: *mut c_char,
    pub mask_bed: *mut c_char,
    pub mask_bed_digest: u32,
    pub host_version: *mut c_char,
    pub cmd_line: *mut c_char,
    pub override_check: c_int,
    pub test_only: c_int,
    pub show_int_params: c_int,
    pub read_buf: *mut u8,
    pub used_read_buf: c_int,
    pub alt_contig_validate: c_int,
    pub auto_detect_validate: c_int,

    pub auto_detect_dir: *const c_char,
    pub pop_alt_contigs_fname: *mut c_char,
    pub pop_alt_liftover_fname: *mut c_char,
    pub pop_snps_input: *mut c_char,
    pub pop_snps_output: *mut c_char,
}

// External C functions
extern "C" {
    /// Set default hash table parameters
    pub fn setDefaultHashParams(
        config: *mut HashTableConfig,
        dir: *const c_char,
        hash_table_type: HashTableType,
    );

    /// Generate hash table from configuration
    /// Returns NULL on success, or error message string on failure
    pub fn generateHashTable(
        config: *mut HashTableConfig,
        argc: c_int,
        argv: *mut *mut c_char,
    ) -> *mut c_char;

    /// Free hash table parameters
    pub fn freeHashParams(config: *mut HashTableConfig);

    /// Standard C library strdup
    pub fn strdup(s: *const c_char) -> *mut c_char;
}

/// Safe wrapper for setDefaultHashParams
///
/// # Safety
/// - `config` must be a valid mutable reference
/// - `dir` must be a valid null-terminated C string pointer
pub unsafe fn set_default_hash_params(
    config: &mut HashTableConfig,
    dir: *const c_char,
    hash_table_type: HashTableType,
) {
    setDefaultHashParams(std::ptr::from_mut(config), dir, hash_table_type);
}

/// Safe wrapper for generateHashTable
///
/// # Safety
/// - `config` must be a valid mutable reference to initialized HashTableConfig
/// - `argv` must point to `argc` valid null-terminated C string pointers
pub unsafe fn generate_hash_table(
    config: &mut HashTableConfig,
    argc: c_int,
    argv: *mut *mut c_char,
) -> *mut c_char {
    generateHashTable(std::ptr::from_mut(config), argc, argv)
}

/// Safe wrapper for freeHashParams
///
/// # Safety
/// - `config` must be a valid mutable reference to previously initialized HashTableConfig
pub unsafe fn free_hash_params(config: &mut HashTableConfig) {
    freeHashParams(std::ptr::from_mut(config));
}
