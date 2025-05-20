use std::error::Error;
use std::io::{stdin, stdout, BufReader, BufWriter, Read, Write};
use std::time::Instant;

use clap::Parser;
use env_logger::Env;
use log::{debug, info};
use rayon::ThreadPoolBuilder;

mod processing;
mod types;

use processing::{
    compute_distances, read_and_parse_fasta, read_and_parse_tabular,
    read_and_parse_tabular_distances, write_distances_to_file, OutputFormat, OutputMode,
};
use types::InputFormat;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    /// The input file or '-' for stdin.
    input: String,

    /// The output file or '-' for stdout.
    output: String,

    /// The format of the input file.
    #[arg(value_enum, short = 'i', long, default_value = "fasta")]
    input_format: InputFormat,

    /// The format of the output file.
    #[arg(value_enum, short = 'o', long, default_value = "tabular")]
    output_format: OutputFormat,

    /// A file with precomputed distances that don't have to be calculated again. The file should be in tabular long format and have the separator as specified by the output-sep flag.
    #[arg(long)]
    precomputed_distances: Option<String>,

    /// The separator character for the input file. Relevant for tabular input files.
    #[arg(long, default_value = "\t")]
    input_sep: char,

    /// The separator character for the output file.
    #[arg(long, default_value = "\t")]
    output_sep: char,

    /// The output mode.
    #[arg(value_enum, short = 'm', long, default_value = "lower-triangle")]
    output_mode: OutputMode,

    /// If set, distance calculations will be stopped when this distance is reached. Useful for large datasets.
    #[arg(short = 'd', long, default_value = None)]
    maxdist: Option<usize>,

    /// Number of threads to use. If not set, all available threads will be used.
    #[arg(short = 't', long, default_value = None)]
    threads: Option<usize>,

    /// Skip the header line of the input file. Relevant for tabular input files.
    #[arg(short = 's', long)]
    skip_header: bool,

    /// Enable verbose mode. Outputs debug messages and calculation times.
    #[arg(short = 'v', long)]
    verbose: bool,
}

fn main() -> Result<(), Box<dyn Error>> {
    let opts: Cli = Cli::parse();
    if opts.verbose {
        env_logger::Builder::from_env(Env::default().default_filter_or("debug")).init();
    } else {
        env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();
    }

    let reader: Box<dyn Read> = if opts.input == "-" {
        Box::new(stdin())
    } else {
        Box::new(std::fs::File::open(&opts.input)?)
    };
    let reader = BufReader::new(reader);

    // print version info
    info!("Version: {}", env!("CARGO_PKG_VERSION"));

    // print command line arguments
    debug!("Cli options: {:?}", opts);

    // Set threads
    match opts.threads {
        Some(threads) => {
            info!("Using {} threads", threads);
            ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .unwrap();
        }
        None => {
            let n = rayon::current_num_threads();
            info!("Using all {n} available threads");
        }
    }

    let start = Instant::now();

    let data_map = match opts.input_format {
        InputFormat::Fasta | InputFormat::FastaAll => {
            read_and_parse_fasta(reader, opts.input_format)?
        }
        InputFormat::Cgmlst | InputFormat::CgmlstHash => {
            read_and_parse_tabular(reader, opts.input_format, opts.input_sep, opts.skip_header)?
        }
    };

    let precomputed_distances = opts
        .precomputed_distances
        .map(|file_name| {
            let file = std::fs::File::open(file_name)?;
            let reader = BufReader::new(file);
            read_and_parse_tabular_distances(reader, opts.output_sep)
        })
        .transpose()?;

    debug!("Reading time: {:?}", start.elapsed());
    let start = Instant::now();

    info!("Computing distances and writing to file: {}", &opts.output);
    let distances = compute_distances(
        &data_map,
        opts.maxdist,
        opts.output_mode,
        precomputed_distances.as_ref(),
    );

    let writer: Box<dyn Write> = if opts.output == "-" {
        Box::new(stdout())
    } else {
        Box::new(std::fs::File::create(&opts.output)?)
    };

    let mut writer = BufWriter::new(writer);

    write_distances_to_file(
        distances,
        &mut writer,
        opts.output_sep,
        opts.output_format,
        data_map.len(),
    )?;

    debug!("Computing + Writing time: {:?}", start.elapsed());
    match opts.maxdist {
        Some(maxdist) => info!("Computed distances with a maximum distance of {}", maxdist),
        None => info!("Computed all distances"),
    }

    info!("Done");

    Ok(())
}
