Build with 
```cargo build -r```

Run with 
```./release/distle --args```


```

Usage: distle [OPTIONS] <INPUT> <OUTPUT> [PRECOMPUTED_DISTANCES]

Arguments:
  <INPUT>
          The input file or '-' for stdin

  <OUTPUT>
          The output file or '-' for stdout

  [PRECOMPUTED_DISTANCES]
          A tabular file with precomputed distances that don't have to be calculated again.

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

  -s, --skip-header
          Skip the header line of the input file. Relevant for tabular input files

  -v, --verbose
          Enable verbose mode. Outputs debug messages and calculation times

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
```
