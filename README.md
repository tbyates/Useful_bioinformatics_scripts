# Useful bioinformatics scripts

#### refer to the doc string in each script for additional information about input/output file formats

## trim_aln.py

#### This script trims a multifasta alignment to the boundaries of the first sequence in the fasta file

### Dependencies

- Biopython

### Usage

`python trim_aln.py -i alignment.fasta -o trimmed_alignment.fasta`


### Full usage
```text
usage: trim_aln.py [-h] -i INPUT -o OUTPUT

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input alignment fasta file
  -o OUTPUT, --output OUTPUT
                        output trimmed alignment file
```

## fisher_chromosome.py

#### This script uses Fisher's exact test to determine if a particular chromosome is overrepresented with genes from a gene set of interest

### Dependencies

- pandas
- scipy

### Usage

`python fisher_chromosome.py -i genes_of_interest.txt -b background_genes.txt -o output_file.txt`


### Full usage
```text
usage: fisher_chromosome.py [-h] -i INPUT -b BACKGROUND -o OUTPUT

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input target genes
  -b BACKGROUND, --background BACKGROUND
                        input background file
  -o OUTPUT, --output OUTPUT
                        output file name
```

## profile_genes_gwas_output.py

#### This script uses gwas output, a user defined window (bp), and gene annotations to profile genes around candidate causal SNPs

### Dependencies

- pandas
- bedtools 
- pybedtools

### Usage

`python profile_genes_gwas_output.py -i gwas_output.txt -p pvalue threshold -a annotation_file.txt -b bed_file_genes.txt -g genome_index.genome -w window size (bp) -o output_file.txt`

### Full usage
```text
usage: profile_genes_gwas_output.py [-h] -i INPUT -p PVALUE -a ANNOTATION_FILE
                                    -b BED_FILE -g GENOME_INDEX -w WINDOW_SIZE
                                    -o OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input gwas output file
  -p PVALUE, --pvalue PVALUE
                        pvalue threshold
  -a ANNOTATION_FILE, --annotation_file ANNOTATION_FILE
                        gene annotation file
  -b BED_FILE, --bed_file BED_FILE
                        gene bed file
  -g GENOME_INDEX, --genome_index GENOME_INDEX
                        genome index
  -w WINDOW_SIZE, --window_size WINDOW_SIZE
                        window size (bp)
  -o OUTPUT, --output OUTPUT
                        output tab delimited file name
```

## at_transversion.py

#### This script reads in a multiple sequence alignment of two sequences in fasta format and calculates the AT transversion ratio

### Dependencies

- Biopython
- pandas

### Usage

`python at_transversion.py -i input_alignment.fasta -o output_at_transversion.txt`

### Full usage

```text
usage: at_transversion.py [-h] -i INPUT -o OUTPUT

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input alignment file
  -o OUTPUT, --output OUTPUT
                        output tab-delimited AT transversion file
```