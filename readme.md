# cgMLST definition from PIRATE pangenomes

```
usage: rackham.py [-h] [--carriage-threshold CARRIAGE_THRESHOLD] [--length-tolerance LENGTH_TOLERANCE] [--max-gaps MAX_GAPS] [--mode {mlst,fsac}] --output OUTPUT pangenome

positional arguments:
  pangenome             Pangenome directory created by PIRATE

optional arguments:
  -h, --help            show this help message and exit
  --carriage-threshold CARRIAGE_THRESHOLD, -t CARRIAGE_THRESHOLD
                        Proportion of genomes a locus must be present in to be included [1.0]
  --length-tolerance LENGTH_TOLERANCE, -l LENGTH_TOLERANCE
                        Maximum allowed value for (min_length / max_length) for alleles of a given locus [0.00]
  --max-gaps MAX_GAPS, -g MAX_GAPS
                        Maximum allowable gaps in a locus [0]
  --mode {mlst,fsac}    Allele header format to use. mlst = '>aspA_1', fsac = '>1' [fsac]
  --output OUTPUT, -o OUTPUT
                        Output directory containing cgMLST alleles

```