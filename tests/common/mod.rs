use std::path::PathBuf;

pub struct TinyFixture {
    pub ref_dir: PathBuf,
    pub fastq_r1: PathBuf,
    pub fastq_r2: PathBuf,
    pub fastq_single: PathBuf,
    pub fastq_interleaved: PathBuf,
}

impl TinyFixture {
    pub fn new() -> Self {
        let base = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("data")
            .join("tiny");
        Self {
            ref_dir: base.join("tiny-2x1Xrepeats.v8"),
            fastq_r1: base.join("1read-1X.fastq"),
            fastq_r2: base.join("1read-1X-r2.fastq"),
            fastq_single: base.join("1read.fastq"),
            fastq_interleaved: base.join("repeat-1X-interleaved.fastq"),
        }
    }

    pub fn assert_present(&self) {
        assert!(self.ref_dir.exists(), "missing tiny reference data");
        assert!(self.fastq_r1.exists(), "missing tiny fastq r1");
        assert!(self.fastq_r2.exists(), "missing tiny fastq r2");
        assert!(self.fastq_single.exists(), "missing tiny fastq single");
        assert!(
            self.fastq_interleaved.exists(),
            "missing tiny fastq interleaved"
        );
    }
}
