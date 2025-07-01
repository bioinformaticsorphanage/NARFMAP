use std::io::Write; fn main() { let mut file = std::fs::File::create("temp_deps/src/main.rs").unwrap(); file.write_all(b"fn main() { println!(\"bitnuc module structure:\"); }").unwrap(); }
