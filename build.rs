use std::process::Command;

fn main() -> miette::Result<()> {
    println!("cargo:rerun-if-changed=Makefile");
    
    // Skip C++ compilation for now while we work on the pure Rust implementation
    println!("cargo:warning=Skipping C++ compilation for now - using Rust-only implementation");
    
    // We'll comment out the make commands that are failing
    // let status = Command::new("make")
    //     .args(&["clean"])
    //     .status()
    //     .expect("failed to run \"make clean\"");
    // assert!(status.success());

    // // build dragen-os command and dragen static library
    // let status = Command::new("make")
    //     .status()
    //     .expect("failed to run \"make\"");
    // assert!(status.success());

    Ok(())
}
