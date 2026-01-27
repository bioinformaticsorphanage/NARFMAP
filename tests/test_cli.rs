mod common;

use common::TinyFixture;
use std::process::Command;
use tempfile::tempdir;

#[test]
#[ignore = "not yet implemented"]
fn align_cli_emits_sam_for_paired_reads() {
    let fixture = TinyFixture::new();
    fixture.assert_present();

    let temp = tempdir().unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_narfmap"))
        .args([
            "align",
            "--ref-dir",
            fixture.ref_dir.to_str().unwrap(),
            "--fastq1",
            fixture.fastq_r1.to_str().unwrap(),
            "--fastq2",
            fixture.fastq_r2.to_str().unwrap(),
            "--output-directory",
            temp.path().to_str().unwrap(),
            "--output-file-prefix",
            "out",
            "--RGID",
            "rg1",
            "--RGSM",
            "sample",
        ])
        .status()
        .unwrap();

    assert!(status.success());

    let sam_path = temp.path().join("out.sam");
    let sam = std::fs::read_to_string(sam_path).unwrap();
    assert!(sam.lines().any(|line| line.starts_with("@SQ")));
}

#[test]
#[ignore = "not yet implemented"]
fn align_cli_handles_interleaved_fastq() {
    let fixture = TinyFixture::new();
    fixture.assert_present();

    let temp = tempdir().unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_narfmap"))
        .args([
            "align",
            "--ref-dir",
            fixture.ref_dir.to_str().unwrap(),
            "--fastq1",
            fixture.fastq_interleaved.to_str().unwrap(),
            "--interleaved",
            "--output-directory",
            temp.path().to_str().unwrap(),
            "--output-file-prefix",
            "out",
        ])
        .status()
        .unwrap();

    assert!(status.success());

    let sam_path = temp.path().join("out.sam");
    let sam = std::fs::read_to_string(sam_path).unwrap();
    assert!(sam.lines().any(|line| line.starts_with("@RG")));
}
