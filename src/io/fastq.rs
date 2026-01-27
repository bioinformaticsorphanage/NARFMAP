//! FASTQ file parsing using noodles

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
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
        self.inner.read_line(&mut seq_line)?;
        let sequence = seq_line.trim().as_bytes().to_vec();

        // Skip + line
        let mut plus_line = String::new();
        self.inner.read_line(&mut plus_line)?;

        // Read quality
        let mut qual_line = String::new();
        self.inner.read_line(&mut qual_line)?;
        let quality = qual_line.trim().as_bytes().to_vec();

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
