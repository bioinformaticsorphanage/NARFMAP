//! SAM file output

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};

use crate::reference::ReferenceSequence;
use crate::Alignment;

/// SAM file writer
pub struct SamWriter {
    inner: BufWriter<File>,
}

impl SamWriter {
    /// Create a new SAM writer
    pub fn create(path: &Path) -> Result<Self> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create SAM file: {}", path.display()))?;
        Ok(Self {
            inner: BufWriter::new(file),
        })
    }

    /// Write SAM header
    pub fn write_header(
        &mut self,
        reference: &ReferenceSequence,
        rgid: &str,
        rgsm: &str,
    ) -> Result<()> {
        // @HD header
        writeln!(self.inner, "@HD\tVN:1.6\tSO:unsorted")?;

        // @SQ sequence dictionary
        for seq in &reference.sequences {
            writeln!(self.inner, "@SQ\tSN:{}\tLN:{}", seq.name, seq.length)?;
        }

        // @RG read group
        writeln!(self.inner, "@RG\tID:{}\tSM:{}", rgid, rgsm)?;

        // @PG program
        writeln!(
            self.inner,
            "@PG\tID:narfmap\tPN:narfmap\tVN:{}",
            env!("CARGO_PKG_VERSION")
        )?;

        Ok(())
    }

    /// Write an alignment record
    pub fn write_alignment(&mut self, aln: &Alignment, rgid: &str) -> Result<()> {
        // Convert sequence to string
        let seq_str = String::from_utf8_lossy(&aln.sequence);

        // Convert quality to string (Phred+33)
        let qual_str: String = aln
            .quality
            .iter()
            .map(|&q| (q.saturating_add(33)) as char)
            .collect();

        writeln!(
            self.inner,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tRG:Z:{}",
            aln.read_name,
            aln.flag,
            aln.ref_name,
            aln.position,
            aln.mapq,
            aln.cigar,
            aln.mate_ref_name,
            aln.mate_position,
            aln.template_length,
            seq_str,
            qual_str,
            rgid
        )?;

        Ok(())
    }

    /// Flush and close the writer
    pub fn finish(mut self) -> Result<()> {
        self.inner.flush()?;
        Ok(())
    }
}
