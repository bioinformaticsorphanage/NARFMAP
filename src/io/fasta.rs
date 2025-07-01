use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use anyhow::{Context, Result};
use seq_io::fasta::{Reader, Record};

use super::sequence::Sequence;

/// FASTA reader that can iterate through sequences in a FASTA file
pub struct FastaReader {
    reader: Reader<BufReader<File>>,
}

impl FastaReader {
    /// Open a FASTA file for reading
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(&path)
            .with_context(|| format!("Failed to open FASTA file: {}", path.as_ref().display()))?;
        
        let reader = Reader::new(BufReader::new(file));
        
        Ok(Self { reader })
    }
    
    /// Read the next sequence from the FASTA file
    pub fn next_sequence(&mut self) -> Result<Option<Sequence>> {
        match self.reader.next() {
            Some(Ok(record)) => {
                let id = String::from_utf8_lossy(&record.head().to_vec()).into_owned();
                let seq_str = String::from_utf8_lossy(&record.seq().to_vec()).into_owned();
                
                let sequence = Sequence::new(id, &seq_str, None)
                    .with_context(|| "Failed to parse FASTA sequence")?;
                
                Ok(Some(sequence))
            },
            Some(Err(e)) => Err(anyhow::Error::new(e).context("Error reading FASTA record")),
            None => Ok(None), // End of file
        }
    }
    
    /// Iterate through all sequences in the FASTA file
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

/// Load a reference genome from a FASTA file
pub fn load_reference<P: AsRef<Path>>(path: P) -> Result<Vec<Sequence>> {
    let mut reader = FastaReader::from_path(path)?;
    reader.iter_sequences().collect()
} 