use std::path::Path;
use anyhow::Result;
use serde::{Deserialize, Serialize};

/// Metadata about reference sequences used in hash table construction
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceMetadata {
    /// List of sequence information
    pub sequences: Vec<SequenceInfo>,
    /// Total reference length (including padding)
    pub total_length: usize,
    /// Total raw length (without padding)
    pub raw_length: usize,
    /// Total non-N bases
    pub non_n_length: usize,
}

/// Information about a single reference sequence
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SequenceInfo {
    /// Sequence index (0-based)
    pub index: usize,
    /// Sequence name/ID from FASTA header
    pub name: String,
    /// Length of the sequence
    pub length: usize,
    /// Start position in the packed reference
    pub start_position: usize,
    /// Number of leading bases trimmed
    pub begin_trim: usize,
    /// Number of trailing bases trimmed
    pub end_trim: usize,
}

impl ReferenceMetadata {
    /// Create new reference metadata from sequences
    pub fn from_sequences(sequences: &[crate::io::sequence::Sequence]) -> Self {
        let mut seq_infos = Vec::new();
        let mut current_position = 0;
        let mut total_raw_length = 0;
        let mut total_non_n = 0;
        
        for (idx, seq) in sequences.iter().enumerate() {
            let seq_len = seq.len();
            total_raw_length += seq_len;
            
            // Count non-N bases
            let non_n_count = seq.seq.to_string()
                .chars()
                .filter(|&c| c != 'N' && c != 'n')
                .count();
            total_non_n += non_n_count;
            
            seq_infos.push(SequenceInfo {
                index: idx,
                name: seq.id.clone(),
                length: seq_len,
                start_position: current_position,
                begin_trim: 0,
                end_trim: 0,
            });
            
            // Add padding between sequences (DRAGMAP uses 1024-byte alignment)
            let padding = if idx < sequences.len() - 1 {
                1024 - (seq_len % 1024)
            } else {
                0
            };
            current_position += seq_len + padding;
        }
        
        Self {
            sequences: seq_infos,
            total_length: current_position,
            raw_length: total_raw_length,
            non_n_length: total_non_n,
        }
    }
    
    /// Save metadata to JSON file
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let file = std::fs::File::create(path)?;
        serde_json::to_writer_pretty(file, self)?;
        Ok(())
    }
    
    /// Load metadata from JSON file
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = std::fs::File::open(path)?;
        let metadata = serde_json::from_reader(file)?;
        Ok(metadata)
    }
    
    /// Write DRAGMAP-compatible configuration lines
    pub fn write_dragmap_config<W: std::io::Write>(&self, writer: &mut W) -> Result<()> {
        writeln!(writer, "reference_sequences  = {}", self.sequences.len())?;
        writeln!(writer, "reference_len        = {}", self.total_length)?;
        writeln!(writer, "reference_len_raw    = {}", self.raw_length)?;
        writeln!(writer, "reference_len_not_n  = {}", self.non_n_length)?;
        
        // Write individual sequence information
        for seq in &self.sequences {
            writeln!(writer, "reference_sequence{}     = '{}'", seq.index, seq.name)?;
            writeln!(writer, "reference_start{}        = {}", seq.index, seq.start_position)?;
            writeln!(writer, "reference_beg_trim{}     = {}", seq.index, seq.begin_trim)?;
            writeln!(writer, "reference_end_trim{}     = {}", seq.index, seq.end_trim)?;
            writeln!(writer, "reference_len{}          = {}", seq.index, seq.length)?;
        }
        
        Ok(())
    }
}