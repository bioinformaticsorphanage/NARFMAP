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
        writeln!(self.inner, "@RG\tID:{rgid}\tSM:{rgsm}")?;

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::ReferenceSequence;
    use tempfile::NamedTempFile;

    fn tiny_ref_dir() -> std::path::PathBuf {
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("data")
            .join("tiny")
            .join("tiny-2x1Xrepeats.v8")
    }

    fn sample_alignment() -> Alignment {
        Alignment {
            read_name: "read1".to_string(),
            flag: 0,
            ref_name: "ref".to_string(),
            position: 1,
            mapq: 60,
            cigar: "4M".to_string(),
            mate_ref_name: "*".to_string(),
            mate_position: 0,
            template_length: 0,
            sequence: b"ACGT".to_vec(),
            quality: vec![0, 1, 2, 3],
            score: 0,
        }
    }

    #[test]
    fn writes_header_and_alignment() {
        let ref_dir = tiny_ref_dir();
        assert!(ref_dir.exists(), "missing tiny reference data");

        let reference = ReferenceSequence::load(&ref_dir).unwrap();
        let temp = NamedTempFile::new().unwrap();
        let path = temp.path().to_path_buf();

        let mut writer = SamWriter::create(&path).unwrap();
        writer.write_header(&reference, "rg1", "sample").unwrap();
        writer.write_alignment(&sample_alignment(), "rg1").unwrap();
        writer.finish().unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        let mut lines = content.lines();

        assert_eq!(lines.next().unwrap(), "@HD\tVN:1.6\tSO:unsorted");

        let sq_line = lines.next().unwrap();
        let expected_len = reference.total_length();
        assert_eq!(sq_line, format!("@SQ\tSN:ref\tLN:{expected_len}"));

        assert_eq!(lines.next().unwrap(), "@RG\tID:rg1\tSM:sample");
        assert!(lines
            .next()
            .unwrap()
            .starts_with("@PG\tID:narfmap\tPN:narfmap\tVN:"));

        let record = lines.next().unwrap();
        assert_eq!(
            record,
            "read1\t0\tref\t1\t60\t4M\t*\t0\t0\tACGT\t!\"#$\tRG:Z:rg1"
        );
    }
}
