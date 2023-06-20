# Quickstart

## Installation

Kaiju can be downloaded and compiled [from source](https://github.com/bioinformatics-centre/kaiju#downloading-and-compiling-kaiju),
or easily installed via the [bioconda channel](https://anaconda.org/bioconda/kaiju):

```
conda install -c bioconda kaiju
```

## Obtaining a Kaiju index

Kaiju requires a index file created from a reference database of protein sequences.

You can either [create a such an index locally](https://github.com/bioinformatics-centre/kaiju#creating-the-reference-database-and-index) or
[download a pre-built index](https://bioinformatics-centre.github.io/kaiju/downloads.html).

For example, to download the Kaiju index for the NCBI BLAST nr database, download the index file with
```
wget https://kaiju.binf.ku.dk/database/kaiju_db_nr_2023-05-10.tgz
```
and unpack the tar archive with:
```
tar xzf kaiju_db_nr_2023-05-10.tgz
```
which will give these 3 files:
```
kaiju_db_nr.fmi
nodes.dmp
names.dmp
```
The Kaiju index itself is in the file `kaiju_db_nr.fmi`, containing the Borrows-Wheeler-Transform and the FM-Index of the protein sequences, wereas `nodes.dmp` and `names.dmp`
contain the taxonomic tree and taxon names from the [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy).


## Running Kaiju
To run Kaiju with the downloaded and unpacked files run:
```
kaiju -t nodes.dmp -f kaiju_db_nr.fmi -i sequencing_reads.fastq.gz
```

For paired-end reads use:
```
kaiju -t nodes.dmp -f kaiju_db_nr.fmi -i sequencing_reads_R1.fastq.gz -j sequencing_reads_R2.fastq
```

Note: The reads must be in the same order in both files!

Kaiju can read input files in FASTQ and FASTA format, which may be gzip-compressed.

By default, Kaiju will print the output to the terminal (STDOUT).
The output can also be written to a file using the `-o` option:
```
kaiju -t nodes.dmp -f kaiju_db.fmi -i sequencing_reads.fastq.gz -o kaiju.out
```

Kaiju can use multiple parallel threads, which can be specified with the `-z` option, e.g. for using 25 parallel threads:
```
kaiju -z 25 -t nodes.dmp -f kaiju_db.fmi -i sequencing_reads.fastq.gz -o kaiju.out
```

Multiple samples can be processed at once using [kaiju-multi](https://github.com/bioinformatics-centre/kaiju#kaiju-multi).

Kaiju has two run modes and several command-line parameters that influence the classification accuracy, see the [original paper](http://www.nature.com/ncomms/2016/160413/ncomms11257/full/ncomms11257.html) and the [README](https://github.com/bioinformatics-centre/kaiju#run-modes).

## Output
Kaiju will print one line for each read or read pair.
The default output format contains three columns separated by tabs:

1. either C or U, indicating whether the read is classified or unclassified.
2. name of the read
3. NCBI taxon identifier of the assigned taxon

Using the option `-v` enables the verbose output, which will print [additional columns](https://github.com/bioinformatics-centre/kaiju#output-format).

The included program [kaiju2table](https://github.com/bioinformatics-centre/kaiju#creating-classification-summary)
converts Kaiju's output file(s) into a summary table for a given taxonomic rank and [kaiju2krona](https://github.com/bioinformatics-centre/kaiju#creating-input-file-for-krona)
creates a file for making a [Krona](https://github.com/marbl/Krona/wiki/KronaTools) visualisation.

