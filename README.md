![Conda Version](https://img.shields.io/conda/vn/bioconda/distle)
![Static Badge](https://img.shields.io/badge/install%20with-docker-important.svg?style=flat&logo=docker&link=https%3A%2F%2Fquay.io%2Frepository%2Fbiocontainers%2Fdistle)
Fast distance matrixes for large-scale pathogen comparisons.

# Installation
```conda install bioconda::distle```


# Quick start

```
distle input_fasta.aln output_distances.tsv

distle --skip-header --input-format cgmlst-hash --output-format phylip input_cgmlst.tsv output_distances.phy
```

# Building
Build with 
```cargo build -r```

# Options

Run with 
```./target/release/distle --help```


```
Usage: distle [OPTIONS] <INPUT> <OUTPUT>

Arguments:
  <INPUT>
          The input file or '-' for stdin

  <OUTPUT>
          The output file or '-' for stdout

Options:
  -i, --input-format <INPUT_FORMAT>
          The format of the input file
          
          [default: fasta]

          Possible values:
          - cgmlst:      A cgmlst table with allele numbers. Optimized for ChewBBACA output
          - cgmlst-hash: A cgmlst table with SHA1 hashes of the nucleotide of the alleles
          - fasta:       An alignment of nucleotide sequences in FASTA format
          - fasta-all:   An alignment of nucleotide sequences in FASTA format. Counts all differences and not just [ACTG]

  -o, --output-format <OUTPUT_FORMAT>
          The format of the output file
          
          [default: tabular]

          Possible values:
          - tabular: Output the distances in a tabular long format
          - phylip:  Output the distances in a Phylip format

      --precomputed-distances <PRECOMPUTED_DISTANCES>
          A file with precomputed distances that don't have to be calculated again. The file should be in tabular long format and have the separator as specified by the output-sep flag

      --input-sep <INPUT_SEP>
          The separator character for the input file. Relevant for tabular input files
          
          [default: "\t"]

      --output-sep <OUTPUT_SEP>
          The separator character for the output file
          
          [default: "\t"]

  -m, --output-mode <OUTPUT_MODE>
          The output mode
          
          [default: lower-triangle]

          Possible values:
          - lower-triangle: Only output the lower triangle of the distance matrix since it is diagonally symmetric
          - full:           Output the full distance matrix

  -d, --maxdist <MAXDIST>
          If set, distance calculations will be stopped when this distance is reached. Useful for large datasets

  -t, --threads <THREADS>
          Number of threads to use. If not set, all available threads will be used

  -s, --skip-header
          Skip the header line of the input file. Relevant for tabular input files

  -v, --verbose
          Enable verbose mode. Outputs debug messages and calculation times

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
```
