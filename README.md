Build with 
```cargo build -r```

Run with 
```./release/distle --args```


```shell
Usage: distle [OPTIONS] <INPUT> <OUTPUT>

Arguments:
  <INPUT>   
  <OUTPUT>  

Options:
  -i, --input-format <INPUT_FORMAT>    [default: fasta] [possible values: chewbbaca, chewbbaca-hash, fasta, fasta-all]
  -o, --output-format <OUTPUT_FORMAT>  [default: tabular] [possible values: tabular, phylip]
      --input-sep <INPUT_SEP>          [default: "\t"]
      --output-sep <OUTPUT_SEP>        [default: "\t"]
  -m, --output-mode <OUTPUT_MODE>      [default: lower-triangle] [possible values: lower-triangle, full]
  -d, --maxdist <MAXDIST>              
  -v, --verbose                        
  -h, --help      
```
