use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use anyhow::{Context, Result};
use seq_io::fastq::{Reader, Record};

use super::sequence::Sequence;

/// FASTQ reader that can iterate through sequences in a FASTQ file
pub struct FastqReader {
    reader: Reader<BufReader<File>>,
}

impl FastqReader {
    /// Open a FASTQ file for reading
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(&path)
            .with_context(|| format!("Failed to open FASTQ file: {}", path.as_ref().display()))?;
        
        let reader = Reader::new(BufReader::new(file));
        
        Ok(Self { reader })
    }
    
    /// Read the next sequence from the FASTQ file
    pub fn next_sequence(&mut self) -> Result<Option<Sequence>> {
        match self.reader.next() {
            Some(Ok(record)) => {
                let id = record.id()
                    .with_context(|| "Invalid UTF-8 in FASTQ record ID")?
                    .to_string();
                let seq_str = String::from_utf8_lossy(&record.seq()).into_owned();
                let qual = Some(record.qual().to_vec());
                
                let sequence = Sequence::new(id, &seq_str, qual)
                    .with_context(|| "Failed to parse FASTQ sequence")?;
                
                Ok(Some(sequence))
            },
            Some(Err(e)) => Err(anyhow::Error::new(e).context("Error reading FASTQ record")),
            None => Ok(None), // End of file
        }
    }
    
    /// Iterate through all sequences in the FASTQ file
    pub fn iter_sequences(&mut self) -> impl Iterator<Item = Result<Sequence>> + '_ {
        std::iter::from_fn(move || {
            match self.next_sequence() {
                Ok(Some(seq)) => Some(Ok(seq)),
                Ok(None) => None,
                Err(e) => Some(Err(e)),
            }
        })
    }
}

/// A paired-end FASTQ reader that can iterate through paired sequences
pub struct PairedFastqReader {
    reader1: FastqReader,
    reader2: FastqReader,
}

impl PairedFastqReader {
    /// Open a pair of FASTQ files for reading
    pub fn from_paths<P: AsRef<Path>>(path1: P, path2: P) -> Result<Self> {
        let reader1 = FastqReader::from_path(&path1)?;
        let reader2 = FastqReader::from_path(&path2)?;
        
        Ok(Self { reader1, reader2 })
    }
    
    /// Read the next pair of sequences from the FASTQ files
    pub fn next_pair(&mut self) -> Result<Option<(Sequence, Sequence)>> {
        match (self.reader1.next_sequence()?, self.reader2.next_sequence()?) {
            (Some(seq1), Some(seq2)) => Ok(Some((seq1, seq2))),
            (None, None) => Ok(None),
            _ => Err(anyhow::anyhow!("Uneven number of reads in paired FASTQ files")),
        }
    }
    
    /// Iterate through all paired sequences in the FASTQ files
    pub fn iter_pairs(&mut self) -> impl Iterator<Item = Result<(Sequence, Sequence)>> + '_ {
        std::iter::from_fn(move || {
            match self.next_pair() {
                Ok(Some(pair)) => Some(Ok(pair)),
                Ok(None) => None,
                Err(e) => Some(Err(e)),
            }
        })
    }
} 