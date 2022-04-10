# Kaiju

Authors:  
Peter Menzel <pmenzel@gmail.com>   
Anders Krogh <krogh@binf.ku.dk>   


Kaiju is a program for the taxonomic classification of high-throughput sequencing
reads, e.g., Illumina or Roche/454, from whole-genome sequencing of
metagenomic DNA. Reads are directly assigned to taxa using the NCBI taxonomy and a
reference database of protein sequences from microbial and viral genomes.

The program is described in [Menzel, P. et al. (2016) Fast and sensitive taxonomic classification for metagenomics with Kaiju. *Nat. Commun.* 7:11257](http://www.nature.com/ncomms/2016/160413/ncomms11257/full/ncomms11257.html) (open access).

Kaiju can be installed locally (see below) or used via a [web server](http://kaiju.binf.ku.dk/).

See the release notes for all releases [here](http://kaiju.binf.ku.dk/index.html#releases).

### License

Copyright (c) 2015-2022 Peter Menzel and Anders Krogh

Kaiju is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Kaiju is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
See the file LICENSE for more details.

You should have received a copy of the GNU General Public License
along with the source code.  If not, see <http://www.gnu.org/licenses/>.


## Downloading and compiling Kaiju
Kaiju can be downloaded directly from GitHub either as a
[compressed archive](https://github.com/bioinformatics-centre/kaiju/archive/master.tar.gz) or
using the git command line client:
```
git clone https://github.com/bioinformatics-centre/kaiju.git
```
This will create the directory `kaiju` in the current directory.

Kaiju is written in C/C++11 for Linux. It uses the zlib library for reading gzip-compressed files.
If not already installed, it is necessary to install the zlib development library, e.g. on Ubuntu using:
```
sudo apt install libz-dev
```

For compiling Kaiju and its associated programs, type:
```
cd kaiju/src
make
```

After compilation, Kaiju's executable files are available in the `kaiju/bin` directory.
You can add this directory to your shell's `$PATH` variable or copy the files to a directory in your PATH.

## Creating the reference database and index

Before classification of reads, Kaiju's database index needs to be built from
the reference protein database.  You can either create a local index based on
the currently available data from GenBank, or download one of the indexes used
by the [Kaiju web server](http://kaiju.binf.ku.dk/).

For creating a local index, the program `kaiju-makedb` in the `bin/` directory
will download a source database and the taxonomy files from the NCBI FTP server,
convert them into a protein database and construct Kaiju's index (the
Burrows-Wheeler transform and the FM-index) in one go.

`kaiju-makedb` needs `curl` and `wget` for downloading the reference databases.

The downloaded files are several GB in size.
It is therefore recommended to run the program in a directory with at least 500 GB of free space.

Example usage:
```
mkdir kaijudb
cd kaijudb
kaiju-makedb -s <DB>
```
The table below lists the available source databases.
Use the database name shown in the first column as argument to option `-s` in `kaiju-makedb`.
The last column denotes the required memory for running Kaiju with the
respective database and for creating the database (in brackets).

| Option | Description | Sequences<sup>\*</sup> | RAM in GB (makedb)<sup>\*</sup> |
| --- | --- | --- | --- |
| `refseq` | Completely assembled and annotated reference genomes of Archaea, Bacteria, and viruses from the NCBI RefSeq database. | 98 M | 67 (84) |
| `progenomes` |  Representative set of genomes from the [proGenomes](http://progenomes.embl.de/) database and viruses from the NCBI RefSeq database. | 41 M | 30 (35) |
| `viruses` |  Only viruses from the NCBI RefSeq database. | 0.6 M | 0.4  (0.5) |
| `plasmids` |  Plasmid sequences from the NCBI RefSeq database. | 3.7 M | 2.2 (4) |
| `fungi` |  [Fungi](ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi) sequences from the NCBI RefSeq database. | 4.4 M | 4.2 (6.4) |
| `nr` | Subset of NCBI BLAST _nr_ database containing all proteins belonging to Archaea, Bacteria and Viruses. | 249 M | 148 (259) |
| `nr_euk` | Like option `-s nr` and additionally include proteins from fungi and microbial eukaryotes, see taxon list in `bin/kaiju-taxonlistEuk.tsv`. | 277 M | 168 (293) |
| `mar` | Protein sequences from all [Mar databases](https://mmp.sfb.uit.no/). Subsets can be chosen by `mar_ref` or `mar_db`. | 41 M | 28 (35) |
| `rvdb` | Protein sequences from [RVDB-prot](https://rvdb-prot.pasteur.fr/) | 10.5 M | 17 (195) |

\* as of April 2022. The databases can also be downloaded from the [web server page](http://kaiju.binf.ku.dk/server).

By default, `kaiju-makedb` uses 5 parallel threads for constructing the index, which can
be changed by using the option `-t`. Note that a higher number of threads
increases the memory usage during index construction, while reducing the number
of threads decreases memory usage.

After `kaiju-makedb` is finished, only the files `kaiju_db_*.fmi`, `nodes.dmp`,
and `names.dmp` are needed to run Kaiju.

### Custom database
It is also possible to make a custom database from a collection of protein sequences.
The format needs to be a FASTA file in which the headers are the numeric NCBI taxon identifiers of the protein sequences,
which can optionally be prefixed by another identifier (e.g. a counter) followed by an underscore, for example:
```
>1_1358
MAQQRRGGFKRRKKVDFIAANKIEVVDYKDTELLKRFISERGKILPRRVTGTSAKNQRKVVNAIKRARVMALLPFVAEDQN
>2_44689
MASTQNIVEEVQKMLDTYDTNKDGEITKAEAVEYFKGKKAFNPERSAIYLFQVYDKDNDGKITIKELAGDIDFDKALKEYKEKQAKSKQQEAEVEEDIEAFILRHNKDDNTDITKDELIQGFKETGAKDPEKSANFILTEMDTNKDGTITVKELRVYYQKVQKLLNPDQ
>3_352472
MKTKSSNNIKKIYYISSILVGIYLCWQIIIQIIFLMDNSIAILEAIGMVVFISVYSLAVAINGWILVGRMKKSSKKAQYEDFYKKMILKSKILLSTIIIVIIVVVVQDIVINFILPQNPQPYVYMIISNFIVGIADSFQMIMVIFVMGELSFKNYFKFKRIEKQKNHIVIGGSSLNSLPVSLPTVKSNESNESNTISINSENNNSKVSTDDTINNVM
>4_91061
MTNPFENDNYTYKVLKNEEGQYSLWPAFLDVPIGWNVVHKEASRNDCLQYVENNWEDLNPKSNQVGKKILVGKR
...
```
The taxon identifiers must be contained in the [NCBI taxonomy files](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz) nodes.dmp and names.dmp.
Then, Kaiju's index is created using the programs `kaiju-mkbwt` and `kaiju-mkfmi`. For example, if the database FASTA file is called `proteins.faa`, then run:
```
kaiju-mkbwt -n 5 -a ACDEFGHIKLMNPQRSTVWY -o proteins proteins.faa
kaiju-mkfmi proteins
```
which creates the file proteins.fmi that is used by Kaiju.
Note that the protein sequences may only contain the uppercase characters of the standard 20 amino acids, all other
characters need to be removed.

## Running Kaiju
Kaiju requires at least three arguments:
```
kaiju -t nodes.dmp -f kaiju_db_*.fmi -i inputfile.fastq
```
Replace `kaiju_db_*.fmi` by the actual `.fmi` file depending on the selected database.
For example, when running `kaiju-makedb -s refseq`, the corresponding index file is `refseq/kaiju_db_refseq.fmi`.

For paired-end reads use `-i firstfile.fastq` and `-j secondfile.fastq`.

The reads must be in the same order in both files. Kaiju will strip suffixes
from the read names by deleting all characters after a `/` or space.  The read
names are then compared between the first and second file and an error is
issued if they are not identical.

Kaiju can read input files in FASTQ and FASTA format, which may also be gzip-compressed.

By default, Kaiju will print the output to the terminal (STDOUT).
The output can also be written to a file using the `-o` option:
```
kaiju -t nodes.dmp -f kaiju_db.fmi -i inputfile.fastq -o kaiju.out
```

Kaiju can use multiple parallel threads, which can be specified with the `-z` option, e.g. for using 25 parallel threads:
```
kaiju -z 25 -t nodes.dmp -f kaiju_db.fmi -i inputfile.fastq -o kaiju.out
```

**kaiju-multi**  
While `kaiju` can only process one input, `kaiju-multi` can take a comma-separated list of input files (and optionally output files) for processing multiple samples at once:
```
kaiju-multi -z 25 -t nodes.dmp -f kaiju_db.fmi -i sample1_R1.fastq,sample2_R1.fastq,sample3_R1.fastq -j sample1_R2.fastq,sample2_R2.fastq,sample3_R2.fastq  -o sample1.out,sample2.out,sample3.out
```
These lists must have the same length.
It's also possible to merge all outputs into one file using output redirection:
```
kaiju-multi -z 25 -t nodes.dmp -f kaiju_db.fmi -i sample1_R1.fastq,sample2_R1.fastq,sample3_R1.fastq -j sample1_R2.fastq,sample2_R2.fastq,sample3_R2.fastq > all_samples.out
```

### Run modes
The default run mode is **Greedy** with three allowed mismatches.
The number of allowed mismatches can be changed using option `-e`.

The run mode can be changed to **MEM** using option `-a`:
```
kaiju -t nodes.dmp -f kaiju_db.fmi -i inputfile.fastq -a mem
```
The cutoffs for minimum required match length and match score can be changed using the options `-m` (default: 11) and `-s` (default: 65).

In Greedy mode, the option `-E` can be used to filter matches by E-value, similar to blastp.
For example, a cutoff of 0.05 can be set by:
```
kaiju -t nodes.dmp -f kaiju_db.fmi -i inputfile.fastq -a greedy -e 5 -E 0.05
```
NB: The thresholds for minimum match length and score are still applied.

If the input sequences are already protein sequences, use option `-p` to disable translation of the input.

Option `-x` enables filtering of query sequences containing
low-complexity regions by using the SEG algorithm from the blast+ package.
It is enabled by default and can be disabled by the `-X` option.  SEG filtering
is always recommended in order to avoid false positive taxon assignments that
are caused by spurious matches due to simple repeat patterns or other
sequencing noise.

### Output format
Kaiju will print one line for each read or read pair.
The default output format contains three columns separated by tabs.
Using the option `-v` enables the verbose output, which will print additional columns:

1. either C or U, indicating whether the read is classified or unclassified.
2. name of the read
3. NCBI taxon identifier of the assigned taxon
4. the length or score of the best match used for classification
5. the taxon identifiers of all database sequences with the best match
6. the accession numbers of all database sequences with the best match
7. matching fragment sequence(s)

NB: Since the _nr_ database aggregates multiple genes of identical sequences, only the first accession number
for each sequence in the __nr__ source file is kept in Kaiju's database and therefore also in the output file.

The number of taxon identifiers (column 5) and accession numbers (column 5) is limited to 20 entries each in
order to reduce large outputs produced by highly abundant protein sequences in _nr_, e.g. from HIV.

## Classification accuracy

The accuracy of the classification depends both on the choice of the reference
database and the chosen options when running Kaiju. These choices also affect
the speed and memory usage of Kaiju.

For highest sensitivity, it is recommended to use the _nr_ database (+eukaryotes)
as a reference database because it is the most comprehensive set of protein
sequences. Alternatively, use proGenomes over Refseq for increased sensitivity.

Greedy run mode yields a higher sensitivity compared with MEM mode.

For fastest classification, use MEM mode and multiple parallel threads
(`-z`); and for lowest memory usage use the proGenomes reference
database. The number of parallel threads has only little impact on memory usage.

Further, the choice of the minimum required match length (`-m`) in MEM mode or
match score (`-s`) in Greedy mode governs the trade-off between sensitivity and
precision of the classification. Please refer to the paper for a discussion on
this topic.

## Helper programs
### Creating input file for Krona
The program `kaiju2krona` can be used to convert Kaiju's tab-separated output file
into a tab-separated text file, which can be imported into [Krona](https://github.com/marbl/Krona/wiki/KronaTools). It requires the `nodes.dmp`
and `names.dmp` files from the NCBI taxonomy for mapping the taxon identifiers from Kaiju's
output to the corresponding taxon names.
```
kaiju2krona -t nodes.dmp -n names.dmp -i kaiju.out -o kaiju.out.krona
```
The file `kaiju.out.krona` can then be imported into Krona and converted into an HTML file using
Krona's `ktImportText` program:
```
ktImportText -o kaiju.out.html kaiju.out.krona
```

### Creating classification summary
The program `kaiju2table` converts Kaiju's output file(s) into a
summary table for a given taxonomic rank, e.g., genus. It requires the
`nodes.dmp` and `names.dmp` files for mapping the taxon identifiers from the third column in the
Kaiju output to the corresponding taxon names.

Basic usage:
```
kaiju2table -t nodes.dmp -n names.dmp -r genus -o kaiju_summary.tsv kaiju.out [kaiju2.out, ...]
```
The program can also filter out taxa with low abundances, e.g. for only showing genera that
comprise at least 1 percent of the total reads:
```
kaiju2table -t nodes.dmp -n names.dmp -r genus -m 1.0 -o kaiju_summary.tsv kaiju.out [kaiju2.out, ...]
```
Similarly, option `-c` can be used to specify the threshold by absolute read count.

Option `-u` disables counting unclassified reads towards the total number of reads when calculating percentages.

Option `-p` will print the full taxon path instead of just the taxon name.

Instead of printing the full taxon path, option `-l` can be used to specify the
ranks to be printed by supplying a comma-separated list, for example:
`-l superkingdom,phylum,class,order,family,genus,species`.

### Adding taxa names to output file
The program `kaiju-addTaxonNames` appends the name that corresponds to the taxon id in
Kaiju's output file as an additional last column to the output.
```
kaiju-addTaxonNames -t nodes.dmp -n names.dmp -i kaiju.out -o kaiju.names.out
```
Option `-u` will omit unclassified reads.  
Option `-p` will print the full taxon path instead of just the taxon name.  
Option `-r` will print the path containing only to the specified ranks. For example,
`-r phylum,genus` will append the names of phylum and genus to the end of each line.

### Merging outputs
The program `kaiju-mergeOutputs` can merge two tab-separated output files in the
column format (see above) used by Kaiju and Kraken. Only the first three columns are used.

The files need to be sorted by the read name in the second column, for example by:
```
sort -k2,2 kaiju.out  >kaiju.out.sort
sort -k2,2 kraken.out >kraken.out.sort
```
Then both files can be merged:
```
kaiju-mergeOutputs -i kaiju.out.sort -j kraken.out.sort -o combined.out -v
```
The shell's process substitution can be used for sorting without creating intermediate files:
```
kaiju-mergeOutputs -i <(sort -k2,2 kaiju.out) -j <(sort -k2,2 kraken.out) -o combined.out -v
```
The output file will be in the same column format as the input files (but only
contain the first three columns) and it will have the same length as the input
files (which also have to be of same length).  In the case of conflicting taxon
identifiers for a classified read in both input files, `kaiju-mergeOutputs` will use the identifier found in the
first input file (specified by `-i`).  This behavior can be changed using the
`-c` option, which can take four possible values:
- `1`: use taxon identifier from the first input file (default)
- `2`: use taxon identifier from the second input file
- `lca`: use the least common ancestor of the taxon identifiers from both files.
- `lowest`: use the lowest ranking of the two taxon identifiers if they are within the same lineage. Otherwise use the LCA.

Options `lca` and `lowest` require the path to the file `nodes.dmp` by using the `-t` option.

When the two tab-separated output files contain the classification score in the 4th column (by running `kaiju -v`), then option `-s` can be used to give precedence to the classification result with the higher score.

### KaijuX and KaijuP

The programs `kaijux` and `kaijup` can be used for finding the best matching
database sequence for each query sequence without taxonomic classification,
i.e., they will just print the name of the database sequence. Thus, both
programs do not use the `nodes.dmp` file containing the taxonomy, but only need
the `.fmi` database file. While `kaijux` takes nucleotide sequences as input
and translates them into the six reading frames like standard `kaiju`,
`kaijup` takes protein sequences as input, which are directly searched in the
database.  All other parameters remain the same as in standard `kaiju`. In case
of paired-end reads, both mates are searched independently.

To build an index for a custom database, all sequences need to be in a single
FASTA file and may only contain the 20 letters from the standard protein
alphabet `ACDEFGHIKLMNPQRSTVWY`.

For example, building the index (the Burrows-Wheeler transform and FM-index) from the
file with the protein sequences `proteins.faa` is done in two steps by the
programs `kaiju-mkbwt` and `kaiju-mkfmi`:
```
kaiju-mkbwt -n 5 -a ACDEFGHIKLMNPQRSTVWY -o proteins proteins.faa
kaiju-mkfmi proteins
```
This will create two intermediate files `proteins.bwt` and `proteins.sa`, and finally
the file `proteins.fmi`, which is used by Kaiju.

The option `-n` for `kaiju-mkbwt` specifies the number of parallel threads. The more
threads are used, the higher the memory consumption becomes.  The option `-e`
for `kaiju-mkbwt` specifies the exponent of the suffix array checkpoint distances and
therefore determines the trade-off between the size of the suffix array and the
speed of the search. The default value is 5.

