use std::collections::HashSet;
use std::error::Error;
use std::io::{BufRead, BufWriter, Write};

use bio::io::fasta;
use clap::ValueEnum;
use rayon::prelude::*;

use crate::types::InputFormat;
use crate::types::SupportedType;

type InputMatrix = Vec<(String, Vec<SupportedType>)>;
#[derive(Debug, PartialEq, Clone, Copy, ValueEnum)]
pub enum OutputMode {
    /// Only output the lower triangle of the distance matrix since it is diagonally symmetric
    LowerTriangle,
    /// Output the full distance matrix
    Full,
}

#[derive(Debug, PartialEq, Clone, Copy, ValueEnum)]
pub enum OutputFormat {
    /// Output the distances in a tabular long format
    Tabular,
    /// Output the distances in a Phylip format
    Phylip,
}

pub fn read_and_parse_tabular<R: BufRead>(
    reader: R,
    input_format: InputFormat,
    separator: char,
    skip_header: bool,
) -> Result<InputMatrix, Box<dyn Error>> {
    // let reader = BufReader::new(reader);
    let mut lines = reader.lines();

    if skip_header {
        lines.next();
    }

    let mut data_vec = Vec::new();

    for line in lines {
        let line = line?;
        let mut fields = line.split(separator);
        let id = match fields.next() {
            Some(id) => String::from(id),
            None => return Err("Missing ID field at start of line".into()),
        };
        let row_data = fields
            .map(|s| SupportedType::from_str(s, input_format))
            .collect::<Result<Vec<_>, _>>()?;
        data_vec.push((id, row_data));
    }

    Ok(data_vec)
}

pub fn read_and_parse_fasta<R: BufRead>(
    reader: R,
    input_format: InputFormat,
) -> Result<InputMatrix, Box<dyn Error>> {
    let reader = fasta::Reader::new(reader);
    let mut data_vec = Vec::new();

    for record in reader.records() {
        let record = record?;
        let id = record.id().to_string();
        let row_data = record
            .seq()
            .iter()
            .map(|&u| SupportedType::from_u8(u, input_format))
            .collect::<Result<Vec<_>, _>>()?;
        data_vec.push((id, row_data));
    }

    Ok(data_vec)
}

pub fn remove_identical_columns<X>(data_map: &mut Vec<(X, Vec<SupportedType>)>) -> usize {
    let (longest_vec_index, longest_len) = data_map
        .iter()
        .enumerate()
        .map(|(i, (_, row))| (i, row.len()))
        .max_by_key(|&(_, len)| len)
        .unwrap_or((0, 0));

    let to_remove: HashSet<_> = (0..longest_len)
        .rev()
        .filter(|&i| {
            let longest = data_map[longest_vec_index].1.get(i);
            data_map.iter().all(|(_, row)| match (longest, row.get(i)) {
                (Some(f), Some(r)) => SupportedType::eq_whithout_exeptions(f, r),
                _ => false,
            })
        })
        .collect();

    for row in data_map {
        let mut index = 0;
        row.1.retain(|_| {
            let keep = !to_remove.contains(&index);
            index += 1;
            keep
        });
    }
    to_remove.len()
}

pub fn compute_distances(
    data_map: &[(String, Vec<impl PartialEq + Clone + std::marker::Sync>)],
    maxdist: Option<usize>,
    output_mode: OutputMode,
) -> impl Iterator<Item = (&str, &str, usize)> + Clone {
    let len = data_map.len();

    (0..len)
        .into_par_iter()
        .flat_map(move |i| {
            let max_j = match output_mode {
                OutputMode::LowerTriangle => i,
                OutputMode::Full => len,
            };
            (0..max_j).into_par_iter().map(move |j| {
                let (id1, row1) = &data_map[i];
                let (id2, row2) = &data_map[j];
                let dist = compute_distance_eq(row1, row2, maxdist);
                (id1.as_str(), id2.as_str(), dist)
            })
        })
        .collect::<Vec<_>>()
        .into_iter()
}

fn compute_distance_eq<T: PartialEq>(row1: &[T], row2: &[T], maxdist: Option<usize>) -> usize {
    row1.iter()
        .zip(row2.iter())
        .filter(|(x, y)| x != y)
        .take(maxdist.unwrap_or(usize::MAX))
        .count()
}

pub fn write_distances_to_file<'a, W: Write>(
    distances: impl Iterator<Item = (&'a str, &'a str, usize)> + Clone,
    writer: W,
    output_sep: char,
    output_format: OutputFormat,
    number_of_samples: usize,
) -> Result<(), Box<dyn Error>> {
    let writer = BufWriter::new(writer);

    match output_format {
        OutputFormat::Tabular => write_distances_to_long_format(distances, writer, output_sep),
        OutputFormat::Phylip => {
            write_distances_to_philip(distances, writer, output_sep, number_of_samples)
        }
    }
}

fn write_distances_to_long_format<'a, W: Write>(
    distances: impl Iterator<Item = (&'a str, &'a str, usize)>,
    mut writer: W,
    output_sep: char,
) -> Result<(), Box<dyn Error>> {
    for (id1, id2, dist) in distances {
        writeln!(writer, "{}{}{}{}{}", id1, output_sep, id2, output_sep, dist)?;
    }
    Ok(())
}

fn write_distances_to_philip<'a, W: Write>(
    distances: impl Iterator<Item = (&'a str, &'a str, usize)> + Clone,
    mut writer: W,
    output_sep: char,
    number_of_samples: usize,
) -> Result<(), Box<dyn Error>> {
    write!(writer, "{}", number_of_samples)?;

    let mut first = true;
    let mut prev_id = "";
    for (id, id2, dist) in distances {
        if first && id != id2 {
            writeln!(writer)?;
            write!(writer, "{}", id2)?;
        }
        if id != prev_id {
            writeln!(writer)?;
            write!(writer, "{}", id)?;
            prev_id = id;
        }
        write!(writer, "{}{}", output_sep, dist)?;
        first = false;
    }
    writeln!(writer)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::types::{ChewBBACAinteger, Nucleotide, NucleotideAll, SHA1Hash};
    use std::str::FromStr;

    use super::*;

    #[test]
    fn test_compute_distance_eq() {
        let row1 = vec![1, 2, 3, 4, 5];
        let row2 = vec![1, 2, 3, 4, 5];
        assert_eq!(compute_distance_eq(&row1, &row2, None), 0);

        let row1 = vec![1, 2, 3, 4, 5];
        let row2 = vec![1, 2, 3, 4, 6];
        assert_eq!(compute_distance_eq(&row1, &row2, None), 1);

        let row1 = vec![1, 2, 3, 4, 5];
        let row2 = vec![1, 2, 3, 4, 6];
        assert_eq!(compute_distance_eq(&row1, &row2, Some(1)), 1);
    }

    #[test]
    fn test_compute_distance_eq_for_chewbbaca() {
        let x0 = ChewBBACAinteger::from_str("-").unwrap();
        let x1 = ChewBBACAinteger::from_str("1").unwrap();
        let x2 = ChewBBACAinteger::from_str("2").unwrap();
        let x3 = ChewBBACAinteger::from_str("3").unwrap();
        let plot = ChewBBACAinteger::from_str("plot5").unwrap();
        let x3_new = ChewBBACAinteger::from_str("INF-3").unwrap();

        let row1 = vec![x0, x1, x2, x3, x1];
        let row2 = vec![x0, x1, x1, x2, x1];
        let row3 = vec![x0, x1, x2, x3_new, plot];
        assert_eq!(compute_distance_eq(&row1, &row3, None), 0);
        assert_eq!(compute_distance_eq(&row2, &row3, None), 2);
        assert_eq!(compute_distance_eq(&row1, &row2, None), 2);
    }

    #[test]
    fn test_compute_distance_eq_for_chewbbaca_hash() {
        let x0 = SHA1Hash::from_str("-").unwrap();
        let x1 = SHA1Hash::from_str("6bc8d04609de559621859873ef301f221cf5d991").unwrap();
        let x2 = SHA1Hash::from_str("1e354c3d41dc0d3c403db19f22de23299a33a1c8").unwrap();
        let x3 = SHA1Hash::from_str("beb636132e9cb496f1c1d37ecafdd62ed02060b0").unwrap();
        let row1 = vec![x0, x1, x2, x3, x1];
        let row2 = vec![x0, x1, x1, x2, x1];
        let row3 = vec![x0, x0, x2, x0, x0];
        assert_eq!(compute_distance_eq(&row1, &row2, None), 2);
        assert_eq!(compute_distance_eq(&row1, &row3, None), 0);
        assert_eq!(compute_distance_eq(&row2, &row3, None), 1);
    }

    #[test]
    fn test_compute_distance_eq_for_fasta() {
        let a = Nucleotide::from_str("A").unwrap();
        let c = Nucleotide::from_str("C").unwrap();
        let g = Nucleotide::from_str("G").unwrap();
        let t = Nucleotide::from_str("T").unwrap();
        let n = Nucleotide::from_str("n").unwrap();
        let d = Nucleotide::from_str("-").unwrap();

        let row1 = vec![a, c, g, t, n, d];
        let row2 = vec![a, n, d, d, d, n];
        let row3 = vec![c, c, g, t, n, d];

        assert_eq!(compute_distance_eq(&row1, &row2, None), 0);
        assert_eq!(compute_distance_eq(&row1, &row3, None), 1);
        assert_eq!(compute_distance_eq(&row2, &row3, None), 1);

        // >ref
        // TACCGTG
        // >sampleA
        // CGTTACT
        // >sampleB
        // NNCNGTN
        let ref1 = vec![t, a, c, c, g, t, g];
        let sample_a = vec![c, g, t, t, a, c, t];
        let sample_b = vec![n, n, c, n, g, t, n];

        assert_eq!(compute_distance_eq(&ref1, &sample_a, None), 7);
        assert_eq!(compute_distance_eq(&ref1, &sample_b, None), 0);
        assert_eq!(compute_distance_eq(&sample_a, &sample_b, None), 3);
    }

    #[test]
    fn test_compute_distance_eq_for_fasta_all() {
        let a = NucleotideAll::from_str("A").unwrap();
        let c = NucleotideAll::from_str("C").unwrap();
        let g = NucleotideAll::from_str("G").unwrap();
        let t = NucleotideAll::from_str("T").unwrap();
        let n = NucleotideAll::from_str("n").unwrap();
        let d = NucleotideAll::from_str("-").unwrap();

        let row1 = vec![a, c, g, t, n, d];
        let row2 = vec![a, n, d, d, d, n];
        let row3 = vec![c, c, g, t, n, d];

        assert_eq!(compute_distance_eq(&row1, &row2, None), 5);
        assert_eq!(compute_distance_eq(&row1, &row3, None), 1);
        assert_eq!(compute_distance_eq(&row2, &row3, None), 6);
    }

    #[test]
    fn test_remove_identical_columns() {
        let a = SupportedType::NucleotideAll(NucleotideAll::from_str("A").unwrap());
        let c = SupportedType::NucleotideAll(NucleotideAll::from_str("C").unwrap());
        let g = SupportedType::NucleotideAll(NucleotideAll::from_str("G").unwrap());
        let t = SupportedType::NucleotideAll(NucleotideAll::from_str("T").unwrap());
        let n = SupportedType::NucleotideAll(NucleotideAll::from_str("n").unwrap());
        let d = SupportedType::NucleotideAll(NucleotideAll::from_str("-").unwrap());

        let mut data_map = vec![
            (String::from("id1"), vec![a, c, g, t, n, d]),
            (String::from("id2"), vec![a, n, d, d, n, n]),
            (String::from("id3"), vec![a, c, g, t, n, d]),
        ];
        let data_map_original = data_map.clone();
        let expected_number_removed = 2;
        let number_removed = remove_identical_columns(&mut data_map);

        let expected_data_map = vec![
            (String::from("id1"), vec![c, g, t, d]),
            (String::from("id2"), vec![n, d, d, n]),
            (String::from("id3"), vec![c, g, t, d]),
        ];

        assert_eq!(number_removed, expected_number_removed);
        assert_eq!(data_map, expected_data_map);

        // test that the distances are also the same
        let distances_original: Vec<_> =
            compute_distances(&data_map_original, None, OutputMode::Full).collect();
        let distances: Vec<_> = compute_distances(&data_map, None, OutputMode::Full).collect();

        assert_eq!(distances_original, distances);
    }
}
