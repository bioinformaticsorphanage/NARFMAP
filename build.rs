use std::process::Command;

fn main() -> miette::Result<()> {
    println!("cargo:rerun-if-changed=Makefile");
    let status = Command::new("make")
        .args(&["clean"])
        .status()
        .expect("failed to run \"make clean\"");
    assert!(status.success());

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
