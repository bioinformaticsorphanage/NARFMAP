mod common;

use common::TinyFixture;
use narfmap::alignment::Aligner;
use narfmap::io::FastqReader;
use narfmap::reference::hash_table::HashTable;
use narfmap::reference::ReferenceSequence;

#[test]
#[ignore = "not yet implemented"]
fn aligns_single_read_from_tiny_fastq() {
    let fixture = TinyFixture::new();
    fixture.assert_present();

    let hash_table = HashTable::load(&fixture.ref_dir).unwrap();
    let reference = ReferenceSequence::load(&fixture.ref_dir).unwrap();
    let aligner = Aligner::new(&hash_table, &reference);

    let mut reader = FastqReader::open(&fixture.fastq_r1).unwrap();
    let read = reader.read_record().unwrap().unwrap();
    let alignment = aligner.align(&read);

    assert!(alignment.flag & 4 == 0);
}

#[test]
#[ignore = "not yet implemented"]
fn aligns_interleaved_pairs_from_tiny_fastq() {
    let fixture = TinyFixture::new();
    fixture.assert_present();

    let hash_table = HashTable::load(&fixture.ref_dir).unwrap();
    let reference = ReferenceSequence::load(&fixture.ref_dir).unwrap();
    let aligner = Aligner::new(&hash_table, &reference);

    let mut reader = FastqReader::open(&fixture.fastq_interleaved).unwrap();
    let read1 = reader.read_record().unwrap().unwrap();
    let read2 = reader.read_record().unwrap().unwrap();
    let aln1 = aligner.align(&read1);
    let aln2 = aligner.align(&read2);

    assert!(aln1.flag & 4 == 0);
    assert!(aln2.flag & 4 == 0);
}
