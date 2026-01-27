use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let manifest_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    let manifest_path = PathBuf::from(&manifest_dir);
    let out_dir = env::var("OUT_DIR").unwrap();
    let out_path = PathBuf::from(&out_dir);

    // Paths to thirdparty dragen code
    let dragen_src = manifest_path.join("thirdparty/dragen/src");
    let hash_gen_dir = dragen_src.join("common/hash_generation");
    let infra_dir = dragen_src.join("host/infra");

    // Collect all C source files for hash generation
    let hash_gen_sources: Vec<PathBuf> = [
        "crc_hash.c",
        "gen_hash_table.c",
        "hash_cfg_file.c",
        "hash_table.c",
        "hash_table_compress.c",
        "host_version_hashtable.c",
        "liftover.c",
        "methylation_hash_table.c",
    ]
    .iter()
    .map(|f| hash_gen_dir.join(f))
    .collect();

    // CRC implementation files
    let crc_sources: Vec<PathBuf> = ["crypto/fast_nonvector_crc32c.c", "crypto/crc32_hw.c"]
        .iter()
        .map(|f| infra_dir.join(f))
        .collect();

    // Include paths
    let includes: Vec<PathBuf> = vec![
        hash_gen_dir.clone(),
        infra_dir.join("public"),
        infra_dir.join("crypto"),
        dragen_src.join("common/public"),
        dragen_src.join("host/dragen_api/sampling"),
        dragen_src.join("host/metrics/public"),
    ];

    // Compile each source file to an object file
    let mut objects = Vec::new();

    for src in hash_gen_sources.iter().chain(crc_sources.iter()) {
        let obj_name = src.file_stem().unwrap().to_str().unwrap();
        let obj_path = out_path.join(format!("{}.o", obj_name));

        let mut cmd = Command::new("clang");
        cmd.arg("-c")
            .arg("-o")
            .arg(&obj_path)
            .arg(src)
            .arg("-std=gnu99")
            .arg("-O2")
            .arg("-DLOCAL_BUILD")
            .arg("-DDRAGEN_OS_BUILD")
            .arg("-D_TARGET_PPC_") // Disable x86-specific rdtsc
            .arg("-Wno-unused-function")
            .arg("-Wno-unused-variable")
            .arg("-Wno-sign-compare");

        for inc in &includes {
            cmd.arg(format!("-I{}", inc.display()));
        }

        let status = cmd
            .status()
            .expect(&format!("Failed to compile {}", src.display()));
        if !status.success() {
            panic!("Compilation failed for {}", src.display());
        }

        objects.push(obj_path);
    }

    // Create static library using libtool (macOS) or ar (Linux)
    let lib_path = out_path.join("libdragen_hash_gen.a");

    #[cfg(target_os = "macos")]
    {
        let mut cmd = Command::new("libtool");
        cmd.arg("-static").arg("-o").arg(&lib_path);
        for obj in &objects {
            cmd.arg(obj);
        }
        let status = cmd.status().expect("Failed to run libtool");
        if !status.success() {
            panic!("libtool failed");
        }
    }

    #[cfg(not(target_os = "macos"))]
    {
        let mut cmd = Command::new("ar");
        cmd.arg("rcs").arg(&lib_path);
        for obj in &objects {
            cmd.arg(obj);
        }
        let status = cmd.status().expect("Failed to run ar");
        if !status.success() {
            panic!("ar failed");
        }
    }

    // Tell cargo to link against our library
    println!("cargo:rustc-link-search=native={}", out_dir);
    println!("cargo:rustc-link-lib=static=dragen_hash_gen");

    // Link against zlib (required for compression)
    println!("cargo:rustc-link-lib=z");

    // Rerun if any source files change
    for src in hash_gen_sources.iter().chain(crc_sources.iter()) {
        println!("cargo:rerun-if-changed={}", src.display());
    }
    println!("cargo:rerun-if-changed=build.rs");
}
