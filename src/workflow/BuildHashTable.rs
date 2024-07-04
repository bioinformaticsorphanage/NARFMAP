#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]

include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

use super::*;
use std::path::PathBuf;

use std::ffi::CString;
use std::mem;
use std::path::Path;

pub fn flow(fasta: PathBuf) {
    println!("==================================================================");
    println!("Building hash table from {fasta:?}");
    println!("==================================================================");
    // unsafe {
    //     // println!("flowing fasta: {:?}", timeval(12, 100));
    //     // generateHashTable(&bhtConfig, opts.argc(), const_cast<char**>(opts.argv()));
    //     generateHashTable()
    // }
    //
    let mut hash_table_types: Vec<HashTableType> = Vec::new();
    hash_table_types.push(0);

    for hash_table_type in &hash_table_types {
        let mut bht_config = mem::MaybeUninit::<hashTableConfig_t>::uninit();
        let bht_config_ptr = bht_config.as_mut_ptr();

        // Set parameters shared with build_hash_table that the user cannot change
        let outdir = Path::new(&opts.output_directory)
            .canonicalize()
            .unwrap()
            .to_str()
            .unwrap();
        setDefaultHashParams(
            bht_config_ptr,
            CString::new(outdir).unwrap().as_ptr(),
            *hash_table_type,
        );
        SetBuildHashTableOptions(&opts, bht_config_ptr, *hash_table_type);

        let error_msg = unsafe {
            let opts_argc = opts.argc();
            let opts_argv = opts
                .argv()
                .iter()
                .map(|s| CString::new(*s).unwrap())
                .collect::<Vec<_>>();
            let opts_argv_ptrs = opts_argv.iter().map(|s| s.as_ptr()).collect::<Vec<_>>();
            generateHashTable(
                bht_config_ptr,
                opts_argc,
                opts_argv_ptrs.as_ptr() as *mut *mut i8,
            )
        };

        if !error_msg.is_null() {
            let error_msg_str = unsafe { CString::from_raw(error_msg).into_string().unwrap() };
            panic!("Hash table generation failed: {}", error_msg_str);
        }

        // FIXME src/lib/workflow/GenHashTableWorkflow.cpp:187
        // We don't need to manage the memory?
        // unsafe {
        //     FreeBuildHashTableOptions(bht_config_ptr);
        // }
    }
}
