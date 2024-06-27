use std::path::PathBuf;

use clap::{Parser, Subcommand};

mod workflow;
use workflow::BuildHashTable;
// #[cxx::bridge]
// pub(crate) mod ffi {
//     unsafe extern "C++" {
//         include!("include/workflow/GenHashTableWorkflow.hpp");
//         include!("include/workflow/Input2SamWorkflow.hpp");
//     }
// }

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    /// Turn debugging information on
    #[arg(short, long, action = clap::ArgAction::Count)]
    debug: u8,

    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Index {
        #[arg(short, long, value_name = "FILE")]
        fasta: PathBuf,
    },

    Align {
        #[arg(short, long, value_name = "FILE")]
        fastq1: PathBuf,
        #[arg(short, long, value_name = "FILE")]
        fastq2: Option<PathBuf>,
    },
}

fn main() {
    let cli = Cli::parse();

    // You can see how many times a particular flag or argument occurred
    // Note, only flags can have multiple occurrences
    match cli.debug {
        0 => println!("Debug mode is off"),
        1 => println!("Debug mode is kind of on"),
        2 => println!("Debug mode is on"),
        _ => println!("Don't be crazy"),
    }

    // You can check for the existence of subcommands, and if found use their
    // matches just as you would the top level cmd
    match cli.command {
        Commands::Index { fasta } => {
            let workflow = BuildHashTable::flow(fasta);
        }
        Commands::Align { fastq1, fastq2 } => {}
    }

    // Continued program logic goes here...
}
