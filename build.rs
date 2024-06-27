use std::process::Command;
extern crate bindgen;

use std::env;
use std::path::PathBuf;

fn main() -> miette::Result<()> {
    println!("cargo:rerun-if-changed=Makefile");
    let status = Command::new("make")
        .args(&["clean"])
        .status()
        .expect("failed to run \"make clean\"");
    assert!(status.success());

    let bindings = bindgen::Builder::default()
        .header("wrapper.h")
        .generate()
        .expect("Unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
    // build dragen-os command and dragen static library
    // let status = Command::new("make")
    //     .status()
    //     .expect("failed to run \"make\"");
    // assert!(status.success());

    // Link libraries
    // println!("cargo:rustc-link-lib=static={}", "dragen");
    // println!("cargo:rustc-link-lib=static={}", "dragen-os");

    Ok(())
}
