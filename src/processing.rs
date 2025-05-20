use core::panic;
use std::collections::HashMap;
use std::error::Error;
use std::io::{BufRead, BufWriter, Write};
use std::str::FromStr;

use bio::io::fasta;
use clap::ValueEnum;
use rayon::prelude::*;

use crate::types::{InputFormat, InputMatrix, SupportedTypeVec};

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
    let mut lines = reader.lines();

    if skip_header {
        lines.next();
    }

    let mut data_vec = Vec::new();

    for line in lines {
        let line = line?;
        let mut fields = line.split(separator);
        let id = fields
            .next()
            .ok_or("Missing ID field at the start of the line")?;
        let id = id.to_string();

        let row_data = match input_format {
            InputFormat::Cgmlst => SupportedTypeVec::Cgmlst(parse_fields(fields)?),
            InputFormat::CgmlstHash => SupportedTypeVec::SHA1Hash(parse_fields(fields)?),
            _ => return Err("Input format not implemented".into()),
        };

        data_vec.push((id, row_data));
    }

    Ok(data_vec)
}

fn parse_fields<'a, I, T>(fields: I) -> Result<Vec<T>, Box<dyn Error>>
where
    I: Iterator<Item = &'a str>,
    T: FromStr,
    T::Err: Error + 'static,
{
    fields
        .map(|s| s.parse::<T>().map_err(|e| Box::new(e) as Box<dyn Error>))
        .collect()
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

        let row_data = match input_format {
            InputFormat::Fasta => SupportedTypeVec::Nucleotide(parse_fasta_seq(record.seq())?),
            InputFormat::FastaAll => {
                SupportedTypeVec::NucleotideAll(parse_fasta_seq(record.seq())?)
            }
            _ => return Err("Input format not implemented".into()),
        };

        data_vec.push((id, row_data));
    }

    Ok(data_vec)
}

fn parse_fasta_seq<T: From<u8>>(seq: &[u8]) -> Result<Vec<T>, Box<dyn Error>> {
    Ok(seq.iter().map(|&u| T::from(u)).collect())
}

pub fn read_and_parse_tabular_distances<R: BufRead>(
    reader: R,
    separator: char,
) -> Result<HashMap<(String, String), usize>, Box<dyn Error>> {
    let mut distances = HashMap::new();
    for line in reader.lines() {
        let line = line?;
        let mut fields = line.split(separator);
        let id1: String = fields
            .next()
            .ok_or("Missing ID field at start of line")?
            .into();
        let id2: String = fields
            .next()
            .ok_or("Missing ID field at start of line")?
            .into();
        let dist = fields.next().ok_or("Missing distance field")?.parse()?;
        distances.insert((id1.clone(), id2.clone()), dist);
        distances.insert((id2, id1), dist); // Also insert the reverse in case the input has a different order
    }
    Ok(distances)
}

// Constants for chunk size calculation
const MIN_CHUNK_SIZE: usize = 100;
const ITEMS_PER_CORE: usize = 100;

pub fn compute_distances<'a>(
    data_map: &'a InputMatrix,
    maxdist: Option<usize>,
    output_mode: OutputMode,
    already_computed: Option<&'a HashMap<(String, String), usize>>,
) -> impl Iterator<Item = (&'a str, &'a str, usize)> + Clone {
    let len = data_map.len();
    let num_cores = rayon::current_num_threads();

    // Calculate chunk size based on number of cores
    // This is optimized to use all available cores while keeping memory usage low
    let chunk_size = len
        .div_ceil(num_cores * ITEMS_PER_CORE)
        .max(MIN_CHUNK_SIZE)
        .min(len / 4)
        .max(1);

    (0..len).step_by(chunk_size).flat_map(move |chunk_start| {
        // Process chunk in parallel and collect results per chunk
        // This is to reduce memory usage while still using all available cores
        let chunk_end = (chunk_start + chunk_size).min(len);

        (chunk_start..chunk_end)
            .into_par_iter()
            .flat_map_iter(|i| {
                let max_j = match output_mode {
                    OutputMode::LowerTriangle => i,
                    OutputMode::Full => len,
                };

                (0..max_j).map(move |j| {
                    let (id1, row1) = &data_map[i];
                    let (id2, row2) = &data_map[j];

                    let dist = already_computed
                        .and_then(|distances| {
                            distances.get(&(id1.to_string(), id2.to_string())).copied()
                        })
                        .unwrap_or_else(|| calculate_distance(row1, row2, maxdist));

                    (id1.as_str(), id2.as_str(), dist)
                })
            })
            .collect::<Vec<_>>()
            .into_iter()
    })
}

fn calculate_distance(
    row1: &SupportedTypeVec,
    row2: &SupportedTypeVec,
    maxdist: Option<usize>,
) -> usize {
    match (row1, row2) {
        (SupportedTypeVec::Nucleotide(r1), SupportedTypeVec::Nucleotide(r2)) => {
            compute_distance_eq(r1, r2, maxdist)
        }
        (SupportedTypeVec::NucleotideAll(r1), SupportedTypeVec::NucleotideAll(r2)) => {
            compute_distance_eq(r1, r2, maxdist)
        }
        (SupportedTypeVec::Cgmlst(r1), SupportedTypeVec::Cgmlst(r2)) => {
            compute_distance_eq(r1, r2, maxdist)
        }
        (SupportedTypeVec::SHA1Hash(r1), SupportedTypeVec::SHA1Hash(r2)) => {
            compute_distance_eq(r1, r2, maxdist)
        }
        _ => panic!("Unsupported type"),
    }
}

fn compute_distance_eq<T: PartialEq>(row1: &[T], row2: &[T], maxdist: Option<usize>) -> usize {
    let maxdist = maxdist.unwrap_or(usize::MAX);
    let mut count = 0;

    for (x, y) in row1.iter().zip(row2.iter()) {
        if x != y {
            count += 1;
            if count >= maxdist {
                break;
            }
        }
    }
    count
}

pub fn write_distances_to_file<'a, W: Write>(
    distances: impl Iterator<Item = (&'a str, &'a str, usize)>,
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
    distances: impl Iterator<Item = (&'a str, &'a str, usize)>,
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
    use crate::types::{ChewBBACAinteger, Hash, Nucleotide, NucleotideAll};
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
        let x0 = Hash::from_str("-").unwrap();
        let x1 = Hash::from_str("6bc8d04609de559621859873ef301f221cf5d991").unwrap();
        let x2 = Hash::from_str("1e354c3d41dc0d3c403db19f22de23299a33a1c8").unwrap();
        let x3 = Hash::from_str("beb636132e9cb496f1c1d37ecafdd62ed02060b0").unwrap();
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
}
