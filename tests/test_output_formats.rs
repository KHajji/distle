use std::fs::File;
use std::io::{BufReader, Cursor, Read, Seek, SeekFrom};

use distle::processing::{
    compute_distances, read_and_parse_fasta, read_and_parse_tabular, write_distances_to_file,
    OutputFormat, OutputMode,
};
use distle::types::InputFormat;

#[test]
pub fn test_output_long() {
    let input = BufReader::new(File::open("tests/data/input.fasta").unwrap());
    let mut output = Cursor::new(Vec::new());
    let input_format = InputFormat::FastaAll;
    let output_format = OutputFormat::Tabular;
    let output_sep = '\t';
    let output_mode = OutputMode::LowerTriangle;
    let maxdist = None;

    let data_map = read_and_parse_fasta(input, input_format).unwrap();
    // remove_identical_columns(&mut data_map);
    let distances = compute_distances(&data_map, maxdist, output_mode, None);
    write_distances_to_file(
        distances,
        &mut output,
        output_sep,
        output_format,
        data_map.len(),
    )
    .unwrap();
    let expected = include_bytes!("data/output.tsv").to_vec();
    let mut result = Vec::new();
    output.seek(SeekFrom::Start(0)).unwrap();
    output.read_to_end(&mut result).unwrap();

    assert_eq!(expected, result);
}

#[test]
pub fn test_output_long_all() {
    let input = BufReader::new(File::open("tests/data/input.fasta").unwrap());
    let mut output = Cursor::new(Vec::new());
    let input_format = InputFormat::FastaAll;
    let output_format = OutputFormat::Tabular;
    let output_sep = '\t';
    let output_mode = OutputMode::Full;
    let maxdist = None;

    let data_map = read_and_parse_fasta(input, input_format).unwrap();
    // remove_identical_columns(&mut data_map);
    let distances = compute_distances(&data_map, maxdist, output_mode, None);
    write_distances_to_file(
        distances,
        &mut output,
        output_sep,
        output_format,
        data_map.len(),
    )
    .unwrap();
    let expected = include_bytes!("data/output_full.tsv").to_vec();
    let mut result = Vec::new();
    output.seek(SeekFrom::Start(0)).unwrap();
    output.read_to_end(&mut result).unwrap();

    assert_eq!(expected, result);
}

#[test]
pub fn test_output_phylip() {
    let input = BufReader::new(File::open("tests/data/input.fasta").unwrap());
    let mut output = Cursor::new(Vec::new());
    let input_format = InputFormat::FastaAll;
    let output_format = OutputFormat::Phylip;
    let output_sep = '\t';
    let output_mode = OutputMode::LowerTriangle;
    let maxdist = None;

    let data_map = read_and_parse_fasta(input, input_format).unwrap();
    // remove_identical_columns(&mut data_map);
    let distances = compute_distances(&data_map, maxdist, output_mode, None);
    write_distances_to_file(
        distances,
        &mut output,
        output_sep,
        output_format,
        data_map.len(),
    )
    .unwrap();
    let expected = include_bytes!("data/output.phylip").to_vec();
    let mut result = Vec::new();
    output.seek(SeekFrom::Start(0)).unwrap();
    output.read_to_end(&mut result).unwrap();

    assert_eq!(expected, result);
}

#[test]
pub fn test_output_phylip_full() {
    let input = BufReader::new(File::open("tests/data/input.fasta").unwrap());
    let mut output = Cursor::new(Vec::new());
    let input_format = InputFormat::FastaAll;
    let output_format = OutputFormat::Phylip;
    let output_sep = '\t';
    let output_mode = OutputMode::Full;
    let maxdist = None;

    let data_map = read_and_parse_fasta(input, input_format).unwrap();
    // remove_identical_columns(&mut data_map);
    let distances = compute_distances(&data_map, maxdist, output_mode, None);
    write_distances_to_file(
        distances,
        &mut output,
        output_sep,
        output_format,
        data_map.len(),
    )
    .unwrap();
    let expected = include_bytes!("data/output_full.phylip").to_vec();
    let mut result = Vec::new();
    output.seek(SeekFrom::Start(0)).unwrap();
    output.read_to_end(&mut result).unwrap();

    assert_eq!(expected, result);
}

#[test]
pub fn test_input_cgmlst_hash() {
    let input = BufReader::new(File::open("tests/data/cgmlst_hash.tsv").unwrap());
    let mut output = Cursor::new(Vec::new());
    let input_format = InputFormat::CgmlstHash;
    let output_format = OutputFormat::Phylip;
    let input_sep = '\t';
    let output_sep = '\t';
    let output_mode = OutputMode::LowerTriangle;
    let maxdist = None;

    let data_map = read_and_parse_tabular(input, input_format, input_sep, false).unwrap();
    // remove_identical_columns(&mut data_map);
    let distances = compute_distances(&data_map, maxdist, output_mode, None);
    write_distances_to_file(
        distances,
        &mut output,
        output_sep,
        output_format,
        data_map.len(),
    )
    .unwrap();
    let expected = include_bytes!("data/output_cgmlst_hash.phylip").to_vec();
    let mut result = Vec::new();
    output.seek(SeekFrom::Start(0)).unwrap();
    output.read_to_end(&mut result).unwrap();

    assert_eq!(expected, result);
}

#[test]
pub fn test_input_cgmlst_hash_full() {
    let input = BufReader::new(File::open("tests/data/cgmlst_hash.tsv").unwrap());
    let mut output = Cursor::new(Vec::new()); // In-memory string
    let input_format = InputFormat::CgmlstHash;
    let output_format = OutputFormat::Phylip;
    let input_sep = '\t';
    let output_sep = '\t';
    let output_mode = OutputMode::Full;
    let maxdist = None;

    let data_map = read_and_parse_tabular(input, input_format, input_sep, false).unwrap();
    // remove_identical_columns(&mut data_map);
    let distances = compute_distances(&data_map, maxdist, output_mode, None);
    write_distances_to_file(
        distances,
        &mut output,
        output_sep,
        output_format,
        data_map.len(),
    )
    .unwrap();

    // let expected = include_bytes!("data/output_cgmlst_hash_full.phylip").to_vec();

    let mut result = Vec::new();
    output.seek(SeekFrom::Start(0)).unwrap(); // Reset the cursor position to the start
    output.read_to_end(&mut result).unwrap(); // Read the output into a vector

    let expected = include_bytes!("data/output_cgmlst_hash_full.phylip").to_vec();
    assert_eq!(expected, result);
}

#[test]
pub fn test_remove_identical() {
    let input = BufReader::new(File::open("tests/data/test_remove_identical.fasta").unwrap());
    let data_map = read_and_parse_fasta(input, InputFormat::Fasta).unwrap();
    let data_map_with_removed_columns = data_map.clone();
    // let _n_removed = remove_identical_columns(&mut data_map_with_removed_columns);

    // assert_eq!(data_map, data_map_with_removed_columns);

    let dist_original: Vec<_> =
        compute_distances(&data_map, None, OutputMode::LowerTriangle, None).collect();
    let dist_removed: Vec<_> = compute_distances(
        &data_map_with_removed_columns,
        None,
        OutputMode::LowerTriangle,
        None,
    )
    .collect();

    assert_eq!(dist_original, dist_removed);

    // let n_removed = remove_identical_columns(&mut data_map_with_removed_columns);
}
