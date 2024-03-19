use distle::processing::{read_and_parse_fasta_file, write_distances_to_file, read_and_parse_tabular_file, OutputFormat, OutputMode, compute_distances};
use distle::types::InputFormat;

#[test]
pub fn test_output_long() {
    let input = "tests/data/input.fasta";
    let output = "tests/data/tmp.tsv";
    let input_format = InputFormat::FastaAll;
    let output_format = OutputFormat::Tabular;
    let output_sep = '\t';
    let output_mode = OutputMode::LowerTriangle;
    let maxdist = None;

    let data_map = read_and_parse_fasta_file(input, input_format).unwrap();
    let distances = compute_distances(
        &data_map,
        maxdist,
        output_mode,
    );
    write_distances_to_file(
        distances,
        output,
        output_sep,
        output_format,
        data_map.len(),
    ).unwrap();
    let expected = std::fs::read_to_string("tests/data/output.tsv").unwrap();
    let result = std::fs::read_to_string(output).unwrap();
    
    std::fs::remove_file(output).unwrap();
    assert_eq!(expected, result);
}

#[test]
pub fn test_output_long_all() {
    let input = "tests/data/input.fasta";
    let output = "tests/data/tmp_full.tsv";
    let input_format = InputFormat::FastaAll;
    let output_format = OutputFormat::Tabular;
    let output_sep = '\t';
    let output_mode = OutputMode::Full;
    let maxdist = None;

    let data_map = read_and_parse_fasta_file(input, input_format).unwrap();
    let distances = compute_distances(
        &data_map,
        maxdist,
        output_mode,
    );
    write_distances_to_file(
        distances,
        output,
        output_sep,
        output_format,
        data_map.len(),
    ).unwrap();
    let expected = std::fs::read_to_string("tests/data/output_full.tsv").unwrap();
    let result = std::fs::read_to_string(output).unwrap();
    
    std::fs::remove_file(output).unwrap();
    assert_eq!(expected, result);
}

#[test]
pub fn test_output_phylip() {
    let input = "tests/data/input.fasta";
    let output = "tests/data/tmp.phylip";
    let input_format = InputFormat::FastaAll;
    let output_format = OutputFormat::Phylip;
    let output_sep = '\t';
    let output_mode = OutputMode::LowerTriangle;
    let maxdist = None;

    let data_map = read_and_parse_fasta_file(input, input_format).unwrap();
    let distances = compute_distances(
        &data_map,
        maxdist,
        output_mode,
    );
    write_distances_to_file(
        distances,
        output,
        output_sep,
        output_format,
        data_map.len(),
    ).unwrap();
    let expected = std::fs::read_to_string("tests/data/output.phylip").unwrap();
    let result = std::fs::read_to_string(output).unwrap();
    
    std::fs::remove_file(output).unwrap();
    assert_eq!(expected, result);
}

#[test]
pub fn test_output_phylip_full() {
    let input = "tests/data/input.fasta";
    let output = "tests/data/tmp_full.phylip";
    let input_format = InputFormat::FastaAll;
    let output_format = OutputFormat::Phylip;
    let output_sep = '\t';
    let output_mode = OutputMode::Full;
    let maxdist = None;

    let data_map = read_and_parse_fasta_file(input, input_format).unwrap();
    let distances = compute_distances(
        &data_map,
        maxdist,
        output_mode,
    );
    write_distances_to_file(
        distances,
        output,
        output_sep,
        output_format,
        data_map.len(),
    ).unwrap();
    let expected = std::fs::read_to_string("tests/data/output_full.phylip").unwrap();
    let result = std::fs::read_to_string(output).unwrap();
    
    std::fs::remove_file(output).unwrap();
    assert_eq!(expected, result);
}

#[test]
pub fn test_input_cgmlst_hash() {
    let input = "tests/data/cgmlst_hash.tsv";
    let output = "tests/data/tmp.phylip";
    let input_format = InputFormat::ChewbbacaHash;
    let output_format = OutputFormat::Phylip;
    let input_sep = '\t';
    let output_sep = '\t';
    let output_mode = OutputMode::LowerTriangle;
    let maxdist = None;

    let data_map = read_and_parse_tabular_file(input, input_format, input_sep).unwrap();
    let distances = compute_distances(
        &data_map,
        maxdist,
        output_mode,
    );
    write_distances_to_file(
        distances,
        output,
        output_sep,
        output_format,
        data_map.len(),
    ).unwrap();
    let expected = std::fs::read_to_string("tests/data/output_cgmlst_hash.phylip").unwrap();
    let result = std::fs::read_to_string(output).unwrap();
    
    std::fs::remove_file(output).unwrap();
    assert_eq!(expected, result);
}

#[test]
pub fn test_input_cgmlst_hash_full() {
    let input = "tests/data/cgmlst_hash.tsv";
    let output = "tests/data/tmp_full.phylip";
    let input_format = InputFormat::ChewbbacaHash;
    let output_format = OutputFormat::Phylip;
    let input_sep = '\t';
    let output_sep = '\t';
    let output_mode = OutputMode::Full;
    let maxdist = None;

    let data_map = read_and_parse_tabular_file(input, input_format, input_sep).unwrap();
    let distances = compute_distances(
        &data_map,
        maxdist,
        output_mode,
    );
    write_distances_to_file(
        distances,
        output,
        output_sep,
        output_format,
        data_map.len(),
    ).unwrap();
    let expected = std::fs::read_to_string("tests/data/output_cgmlst_hash_full.phylip").unwrap();
    let result = std::fs::read_to_string(output).unwrap();
    
    std::fs::remove_file(output).unwrap();
    assert_eq!(expected, result);
}
