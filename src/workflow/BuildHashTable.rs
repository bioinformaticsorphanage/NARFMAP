#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

use super::*;
use std::path::PathBuf;
pub fn flow(fasta: PathBuf) {
    unsafe {
        println!("flowing fasta: {:?}", timeval(12, 100));
    }
}
