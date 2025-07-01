use std::io::{BufWriter, Write};
use std::fs::File;
use std::path::Path;
use anyhow::{Result, anyhow};
use log::{info, debug};

use super::AlignmentResult;

/// SAM format flags
pub mod sam_flags {
    pub const PAIRED: u16 = 0x1;
    pub const PROPER_PAIR: u16 = 0x2;
    pub const UNMAPPED: u16 = 0x4;
    pub const MATE_UNMAPPED: u16 = 0x8;
    pub const REVERSE: u16 = 0x10;
    pub const MATE_REVERSE: u16 = 0x20;
    pub const FIRST_IN_PAIR: u16 = 0x40;
    pub const SECOND_IN_PAIR: u16 = 0x80;
    pub const SECONDARY: u16 = 0x100;
    pub const QUALITY_FAIL: u16 = 0x200;
    pub const DUPLICATE: u16 = 0x400;
    pub const SUPPLEMENTARY: u16 = 0x800;
}

/// SAM file writer
pub struct SamWriter {
    writer: BufWriter<Box<dyn Write>>,
}

impl SamWriter {
    /// Create a new SAM writer
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::create(path)?;
        let writer = BufWriter::new(Box::new(file) as Box<dyn Write>);
        Ok(Self { writer })
    }

    /// Create a new SAM writer that outputs to stdout
    pub fn new_stdout() -> Result<Self> {
        use std::io::stdout;
        let writer = BufWriter::new(Box::new(stdout()) as Box<dyn Write>);
        Ok(Self { writer })
    }

    /// Write SAM header
    pub fn write_header(&mut self, reference_sequences: &[(String, u32)]) -> Result<()> {
        // Write version header
        writeln!(self.writer, "@HD\tVN:1.6\tSO:unsorted")?;
        
        // Write reference sequence headers
        for (ref_name, ref_len) in reference_sequences {
            writeln!(self.writer, "@SQ\tSN:{}\tLN:{}", ref_name, ref_len)?;
        }
        
        // Write program header
        writeln!(self.writer, "@PG\tID:narfmap\tPN:narfmap\tVN:{}", env!("CARGO_PKG_VERSION"))?;
        
        Ok(())
    }

    /// Write an aligned read
    pub fn write_alignment(&mut self, alignment: &AlignmentResult, read_sequence: &str, read_quality: &str) -> Result<()> {
        let flags = self.calculate_flags(alignment);
        
        let mate_ref = alignment.mate_reference_id.as_deref().unwrap_or("*");
        let mate_pos = alignment.mate_position.unwrap_or(0);
        let tlen = alignment.template_length.unwrap_or(0);
        
        writeln!(
            self.writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            alignment.read_id,          // QNAME
            flags,                      // FLAG
            alignment.reference_id,     // RNAME
            alignment.position + 1,     // POS (1-based)
            alignment.mapq,             // MAPQ
            alignment.cigar,            // CIGAR
            mate_ref,                   // MRNM
            mate_pos + 1,               // MPOS (1-based)
            tlen,                       // TLEN
            read_sequence,              // SEQ
            read_quality                // QUAL
        )?;
        
        Ok(())
    }

    /// Write an unmapped read
    pub fn write_unmapped(&mut self, read_id: &str, read_sequence: &str, read_quality: &str) -> Result<()> {
        writeln!(
            self.writer,
            "{}\t{}\t*\t0\t0\t*\t*\t0\t0\t{}\t{}",
            read_id,                    // QNAME
            sam_flags::UNMAPPED,        // FLAG
            read_sequence,              // SEQ
            read_quality                // QUAL
        )?;
        
        Ok(())
    }

    /// Calculate SAM flags for an alignment
    fn calculate_flags(&self, alignment: &AlignmentResult) -> u16 {
        let mut flags = 0u16;
        
        if alignment.is_paired {
            flags |= sam_flags::PAIRED;
        }
        
        if alignment.is_proper_pair {
            flags |= sam_flags::PROPER_PAIR;
        }
        
        if alignment.is_reverse {
            flags |= sam_flags::REVERSE;
        }
        
        // TODO: Add more flag logic for paired-end reads
        
        flags
    }

    /// Flush the writer
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush()?;
        Ok(())
    }
}

impl Drop for SamWriter {
    fn drop(&mut self) {
        let _ = self.flush();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;
    use tempfile::NamedTempFile;

    #[test]
    fn test_sam_writer() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        let temp_path = temp_file.path().to_path_buf();
        
        // Create writer and write header
        {
            let mut writer = SamWriter::new(&temp_path)?;
            
            let ref_seqs = vec![
                ("chr1".to_string(), 1000),
                ("chr2".to_string(), 2000),
            ];
            writer.write_header(&ref_seqs)?;
            
            // Write an alignment
            let alignment = AlignmentResult {
                read_id: "read1".to_string(),
                reference_id: "chr1".to_string(),
                position: 100,
                cigar: "50M".to_string(),
                mapq: 60,
                is_reverse: false,
                is_paired: false,
                is_proper_pair: false,
                mate_reference_id: None,
                mate_position: None,
                template_length: None,
            };
            
            writer.write_alignment(&alignment, "ACGTACGTACGT", "############")?;
            writer.write_unmapped("unmapped_read", "TTTTTTTTTTTT", "############")?;
        }
        
        // Read back and verify
        let mut content = String::new();
        temp_file.read_to_string(&mut content)?;
        
        assert!(content.contains("@HD\tVN:1.6"));
        assert!(content.contains("@SQ\tSN:chr1\tLN:1000"));
        assert!(content.contains("read1\t0\tchr1\t101\t60\t50M"));
        assert!(content.contains("unmapped_read\t4\t*\t0\t0\t*"));
        
        Ok(())
    }
}