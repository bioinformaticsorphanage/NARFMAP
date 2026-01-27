//! FASTQ file parsing using noodles

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{bail, Context, Result};
use flate2::read::GzDecoder;

use crate::Read;

/// FASTQ file reader supporting gzip compression
pub struct FastqReader {
    inner: Box<dyn BufRead>,
}

impl FastqReader {
    /// Open a FASTQ file (detects gzip by extension)
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open FASTQ: {}", path.display()))?;

        let inner: Box<dyn BufRead> = if path
            .extension()
            .is_some_and(|ext| ext == "gz" || ext == "gzip")
        {
            Box::new(BufReader::new(GzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        Ok(Self { inner })
    }

    /// Read the next record
    pub fn read_record(&mut self) -> Result<Option<Read>> {
        let mut name_line = String::new();
        let bytes_read = self.inner.read_line(&mut name_line)?;
        if bytes_read == 0 {
            return Ok(None);
        }

        // Parse name (skip @ prefix)
        let name = name_line
            .trim()
            .strip_prefix('@')
            .unwrap_or(&name_line)
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();

        // Read sequence
        let mut seq_line = String::new();
        if self.inner.read_line(&mut seq_line)? == 0 {
            bail!("Incomplete FASTQ record: missing sequence line");
        }
        let sequence = seq_line.trim().as_bytes().to_vec();

        // Skip + line
        let mut plus_line = String::new();
        if self.inner.read_line(&mut plus_line)? == 0 {
            bail!("Incomplete FASTQ record: missing '+' line");
        }
        if !plus_line.trim_start().starts_with('+') {
            bail!("Invalid FASTQ record: expected '+' line");
        }

        // Read quality
        let mut qual_line = String::new();
        if self.inner.read_line(&mut qual_line)? == 0 {
            bail!("Incomplete FASTQ record: missing quality line");
        }
        let quality = qual_line.trim().as_bytes().to_vec();

        if sequence.len() != quality.len() {
            bail!("FASTQ sequence/quality length mismatch");
        }

        Ok(Some(Read {
            name,
            sequence,
            quality,
        }))
    }
}

impl Iterator for FastqReader {
    type Item = Result<Read>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.read_record() {
            Ok(Some(read)) => Some(Ok(read)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{BufReader, Cursor};

    fn reader_from_str(data: &str) -> FastqReader {
        let bytes = data.as_bytes().to_vec();
        FastqReader {
            inner: Box::new(BufReader::new(Cursor::new(bytes))),
        }
    }

    #[test]
    fn read_record_parses_valid_record() {
        let data = "@read1 extra\nACGT\n+\n!!!!\n";
        let mut reader = reader_from_str(data);
        let read = reader.read_record().unwrap().unwrap();

        assert_eq!(read.name, "read1");
        assert_eq!(read.sequence, b"ACGT");
        assert_eq!(read.quality, b"!!!!");
        assert!(reader.read_record().unwrap().is_none());
    }

    #[test]
    fn read_record_errors_on_missing_plus_line() {
        let data = "@read1\nACGT\n*\n!!!!\n";
        let mut reader = reader_from_str(data);

        assert!(reader.read_record().is_err());
    }

    #[test]
    fn read_record_errors_on_incomplete_record() {
        let data = "@read1\nACGT\n+\n";
        let mut reader = reader_from_str(data);

        assert!(reader.read_record().is_err());
    }

    #[test]
    fn read_record_errors_on_length_mismatch() {
        let data = "@read1\nACGT\n+\n!!!\n";
        let mut reader = reader_from_str(data);

        assert!(reader.read_record().is_err());
    }
}
