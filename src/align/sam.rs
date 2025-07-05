use std::io::{BufWriter, Write};
use std::fs::File;
use std::path::Path;
use std::collections::HashMap;
use anyhow::Result;
use chrono::{DateTime, Utc};
// use log::{info, debug}; // Unused for now

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

/// SAM header information
#[derive(Debug, Clone)]
pub struct SamHeader {
    /// Version and sorting order
    pub version: String,
    pub sort_order: String,
    pub group_order: Option<String>,
    
    /// Reference sequences
    pub references: Vec<ReferenceInfo>,
    
    /// Read groups
    pub read_groups: Vec<ReadGroup>,
    
    /// Program information
    pub programs: Vec<ProgramInfo>,
    
    /// Comments
    pub comments: Vec<String>,
}

/// Reference sequence information
#[derive(Debug, Clone)]
pub struct ReferenceInfo {
    pub name: String,
    pub length: u32,
    pub assembly: Option<String>,
    pub md5: Option<String>,
    pub species: Option<String>,
    pub uri: Option<String>,
}

/// Read group information
#[derive(Debug, Clone)]
pub struct ReadGroup {
    pub id: String,
    pub sample: String,
    pub library: Option<String>,
    pub platform: Option<String>,
    pub platform_unit: Option<String>,
    pub center: Option<String>,
    pub description: Option<String>,
    pub date: Option<DateTime<Utc>>,
    pub flow_order: Option<String>,
    pub key_sequence: Option<String>,
    pub predicted_median_insert_size: Option<u32>,
}

/// Program information
#[derive(Debug, Clone)]
pub struct ProgramInfo {
    pub id: String,
    pub name: String,
    pub version: String,
    pub command_line: Option<String>,
    pub previous_program: Option<String>,
    pub description: Option<String>,
}

impl Default for SamHeader {
    fn default() -> Self {
        Self {
            version: "1.6".to_string(),
            sort_order: "unsorted".to_string(),
            group_order: None,
            references: Vec::new(),
            read_groups: Vec::new(),
            programs: vec![ProgramInfo {
                id: "narfmap".to_string(),
                name: "narfmap".to_string(),
                version: env!("CARGO_PKG_VERSION").to_string(),
                command_line: None,
                previous_program: None,
                description: Some("NARFMAP - A Rust implementation of the Dragen mapper/aligner".to_string()),
            }],
            comments: Vec::new(),
        }
    }
}

/// SAM file writer
pub struct SamWriter {
    writer: BufWriter<Box<dyn Write>>,
    header: SamHeader,
}

impl SamWriter {
    /// Create a new SAM writer
    pub fn new<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::create(path)?;
        let writer = BufWriter::new(Box::new(file) as Box<dyn Write>);
        let header = SamHeader::default();
        Ok(Self { writer, header })
    }

    /// Create a new SAM writer that outputs to stdout
    pub fn new_stdout() -> Result<Self> {
        use std::io::stdout;
        let writer = BufWriter::new(Box::new(stdout()) as Box<dyn Write>);
        let header = SamHeader::default();
        Ok(Self { writer, header })
    }

    /// Create a new SAM writer with custom header
    pub fn new_with_header<P: AsRef<Path>>(path: P, header: SamHeader) -> Result<Self> {
        let file = File::create(path)?;
        let writer = BufWriter::new(Box::new(file) as Box<dyn Write>);
        Ok(Self { writer, header })
    }

    /// Create a new SAM writer to stdout with custom header
    pub fn new_stdout_with_header(header: SamHeader) -> Result<Self> {
        use std::io::stdout;
        let writer = BufWriter::new(Box::new(stdout()) as Box<dyn Write>);
        Ok(Self { writer, header })
    }

    /// Add a reference sequence to the header
    pub fn add_reference(&mut self, name: String, length: u32) {
        self.header.references.push(ReferenceInfo {
            name,
            length,
            assembly: None,
            md5: None,
            species: None,
            uri: None,
        });
    }

    /// Add a read group to the header
    pub fn add_read_group(&mut self, read_group: ReadGroup) {
        self.header.read_groups.push(read_group);
    }

    /// Set command line for the program
    pub fn set_command_line(&mut self, command_line: String) {
        if let Some(program) = self.header.programs.first_mut() {
            program.command_line = Some(command_line);
        }
    }

    /// Write SAM header (legacy method for compatibility)
    pub fn write_header(&mut self, reference_sequences: &[(String, u32)]) -> Result<()> {
        // Add reference sequences to header if not already present
        for (ref_name, ref_len) in reference_sequences {
            if !self.header.references.iter().any(|r| r.name == *ref_name) {
                self.add_reference(ref_name.clone(), *ref_len);
            }
        }
        
        self.write_complete_header()
    }

    /// Write complete SAM header with all information
    pub fn write_complete_header(&mut self) -> Result<()> {
        // Write @HD header line
        write!(self.writer, "@HD\tVN:{}\tSO:{}", self.header.version, self.header.sort_order)?;
        if let Some(ref group_order) = self.header.group_order {
            write!(self.writer, "\tGO:{}", group_order)?;
        }
        writeln!(self.writer)?;
        
        // Write @SQ header lines for reference sequences
        for reference in &self.header.references {
            write!(self.writer, "@SQ\tSN:{}\tLN:{}", reference.name, reference.length)?;
            if let Some(ref assembly) = reference.assembly {
                write!(self.writer, "\tAS:{}", assembly)?;
            }
            if let Some(ref md5) = reference.md5 {
                write!(self.writer, "\tM5:{}", md5)?;
            }
            if let Some(ref species) = reference.species {
                write!(self.writer, "\tSP:{}", species)?;
            }
            if let Some(ref uri) = reference.uri {
                write!(self.writer, "\tUR:{}", uri)?;
            }
            writeln!(self.writer)?;
        }
        
        // Write @RG header lines for read groups
        for read_group in &self.header.read_groups {
            write!(self.writer, "@RG\tID:{}\tSM:{}", read_group.id, read_group.sample)?;
            if let Some(ref library) = read_group.library {
                write!(self.writer, "\tLB:{}", library)?;
            }
            if let Some(ref platform) = read_group.platform {
                write!(self.writer, "\tPL:{}", platform)?;
            }
            if let Some(ref platform_unit) = read_group.platform_unit {
                write!(self.writer, "\tPU:{}", platform_unit)?;
            }
            if let Some(ref center) = read_group.center {
                write!(self.writer, "\tCN:{}", center)?;
            }
            if let Some(ref description) = read_group.description {
                write!(self.writer, "\tDS:{}", description)?;
            }
            if let Some(ref date) = read_group.date {
                write!(self.writer, "\tDT:{}", date.format("%Y-%m-%d"))?;
            }
            if let Some(ref flow_order) = read_group.flow_order {
                write!(self.writer, "\tFO:{}", flow_order)?;
            }
            if let Some(ref key_sequence) = read_group.key_sequence {
                write!(self.writer, "\tKS:{}", key_sequence)?;
            }
            if let Some(ref predicted_median_insert_size) = read_group.predicted_median_insert_size {
                write!(self.writer, "\tPI:{}", predicted_median_insert_size)?;
            }
            writeln!(self.writer)?;
        }
        
        // Write @PG header lines for programs
        for program in &self.header.programs {
            write!(self.writer, "@PG\tID:{}\tPN:{}\tVN:{}", program.id, program.name, program.version)?;
            if let Some(ref command_line) = program.command_line {
                write!(self.writer, "\tCL:{}", command_line)?;
            }
            if let Some(ref previous_program) = program.previous_program {
                write!(self.writer, "\tPP:{}", previous_program)?;
            }
            if let Some(ref description) = program.description {
                write!(self.writer, "\tDS:{}", description)?;
            }
            writeln!(self.writer)?;
        }
        
        // Write @CO header lines for comments
        for comment in &self.header.comments {
            writeln!(self.writer, "@CO\t{}", comment)?;
        }
        
        Ok(())
    }

    /// Write an aligned read
    pub fn write_alignment(&mut self, alignment: &AlignmentResult, read_sequence: &str, read_quality: &str) -> Result<()> {
        self.write_alignment_with_tags(alignment, read_sequence, read_quality, &HashMap::new())
    }

    /// Write an aligned read with optional tags
    pub fn write_alignment_with_tags(
        &mut self, 
        alignment: &AlignmentResult, 
        read_sequence: &str, 
        read_quality: &str,
        tags: &HashMap<String, String>
    ) -> Result<()> {
        let flags = self.calculate_flags(alignment);
        
        let mate_ref = alignment.mate_reference_id.as_deref().unwrap_or("*");
        let mate_pos = alignment.mate_position.unwrap_or(0);
        let tlen = alignment.template_length.unwrap_or(0);
        
        // Write the 11 mandatory SAM fields
        write!(
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
        
        // Add read group tag if available
        if !self.header.read_groups.is_empty() {
            write!(self.writer, "\tRG:Z:{}", self.header.read_groups[0].id)?;
        }
        
        // Add any additional tags
        for (tag, value) in tags {
            write!(self.writer, "\t{}:{}", tag, value)?;
        }
        
        writeln!(self.writer)?;
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

    #[test]
    fn test_sam_writer_with_complete_header() -> Result<()> {
        let mut temp_file = NamedTempFile::new()?;
        let temp_path = temp_file.path().to_path_buf();
        
        // Create writer with comprehensive header
        {
            let mut header = SamHeader::default();
            header.sort_order = "coordinate".to_string();
            header.group_order = Some("query".to_string());
            header.comments.push("Generated by NARFMAP test suite".to_string());
            
            // Add reference sequences with full metadata
            header.references = vec![
                ReferenceInfo {
                    name: "chr1".to_string(),
                    length: 248956422,
                    assembly: Some("GRCh38".to_string()),
                    md5: Some("6aef897c3d6ff0c78aff06ac189178dd".to_string()),
                    species: Some("Homo sapiens".to_string()),
                    uri: Some("http://www.ncbi.nlm.nih.gov/nuccore/NC_000001.11".to_string()),
                },
                ReferenceInfo {
                    name: "chr2".to_string(),
                    length: 242193529,
                    assembly: Some("GRCh38".to_string()),
                    md5: Some("f98db672eb0993dcfdabafe2a882905c".to_string()),
                    species: Some("Homo sapiens".to_string()),
                    uri: Some("http://www.ncbi.nlm.nih.gov/nuccore/NC_000002.12".to_string()),
                },
            ];
            
            // Add read groups with comprehensive metadata
            header.read_groups = vec![
                ReadGroup {
                    id: "sample1".to_string(),
                    sample: "NA12878".to_string(),
                    library: Some("lib1".to_string()),
                    platform: Some("ILLUMINA".to_string()),
                    platform_unit: Some("flowcell1.lane1".to_string()),
                    center: Some("BI".to_string()),
                    description: Some("Whole genome sequencing".to_string()),
                    date: Some(chrono::Utc::now()),
                    flow_order: Some("TACG".to_string()),
                    key_sequence: Some("TCAG".to_string()),
                    predicted_median_insert_size: Some(500),
                },
            ];
            
            // Update program info
            header.programs[0].command_line = Some("narfmap align -r ref_dir -1 reads.fastq".to_string());
            
            let mut writer = SamWriter::new_with_header(&temp_path, header)?;
            
            // Write the complete header
            writer.write_complete_header()?;
            
            // Write alignments with read group tags
            let alignment = AlignmentResult {
                read_id: "read1".to_string(),
                reference_id: "chr1".to_string(),
                position: 1000000,
                cigar: "100M".to_string(),
                mapq: 60,
                is_reverse: false,
                is_paired: true,
                is_proper_pair: true,
                mate_reference_id: Some("chr1".to_string()),
                mate_position: Some(1000150),
                template_length: Some(250),
            };
            
            let mut tags = HashMap::new();
            tags.insert("NM:i".to_string(), "0".to_string());
            tags.insert("AS:i".to_string(), "100".to_string());
            tags.insert("XS:i".to_string(), "50".to_string());
            
            writer.write_alignment_with_tags(&alignment, "A".repeat(100).as_str(), "#".repeat(100).as_str(), &tags)?;
        }
        
        // Read back and verify comprehensive header
        let mut content = String::new();
        temp_file.read_to_string(&mut content)?;
        
        // Verify header sections
        assert!(content.contains("@HD\tVN:1.6\tSO:coordinate\tGO:query"));
        assert!(content.contains("@SQ\tSN:chr1\tLN:248956422\tAS:GRCh38\tM5:6aef897c3d6ff0c78aff06ac189178dd\tSP:Homo sapiens\tUR:http://www.ncbi.nlm.nih.gov/nuccore/NC_000001.11"));
        assert!(content.contains("@SQ\tSN:chr2\tLN:242193529"));
        assert!(content.contains("@RG\tID:sample1\tSM:NA12878\tLB:lib1\tPL:ILLUMINA\tPU:flowcell1.lane1\tCN:BI\tDS:Whole genome sequencing"));
        assert!(content.contains("@PG\tID:narfmap\tPN:narfmap\tVN:1.4.2\tCL:narfmap align -r ref_dir -1 reads.fastq"));
        assert!(content.contains("@CO\tGenerated by NARFMAP test suite"));
        
        // Verify alignment with read group and tags (flag 3 = paired + proper pair)
        assert!(content.contains("read1\t3\tchr1\t1000001\t60\t100M\tchr1\t1000151\t250"));
        assert!(content.contains("RG:Z:sample1"));
        assert!(content.contains("NM:i:0"));
        assert!(content.contains("AS:i:100"));
        assert!(content.contains("XS:i:50"));
        
        Ok(())
    }
}