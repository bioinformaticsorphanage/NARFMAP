use seq_io::fasta::Record; fn main() { println!("Checking Record.id() signature"); let record_type = std::any::type_name::<fn() -> Result<&str, _>>(); println!("{}", record_type); }
