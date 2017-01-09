#Kaiju

Authors:  
Peter Menzel <pmenzel@gmail.com>   
Anders Krogh <krogh@binf.ku.dk>   


Kaiju is a program for the taxonomic classification of high-throughput sequencing
reads, e.g., Illumina or Roche/454, from whole-genome sequencing of
metagenomic DNA. Reads are directly assigned to taxa using the NCBI taxonomy and a
reference database of protein sequences from microbial and viral genomes.

The program is described in [Menzel, P. et al. (2016) Fast and sensitive taxonomic classification for metagenomics with Kaiju. *Nat. Commun.* 7:11257](http://www.nature.com/ncomms/2016/160413/ncomms11257/full/ncomms11257.html) (open access).

Kaiju can be installed locally (see below) or used via a [web server](http://kaiju.binf.ku.dk/).

###License

Copyright (c) 2015, 2016 Peter Menzel and Anders Krogh

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


##Downloading and compiling Kaiju
Kaiju can be downloaded directly from GitHub either as a
[compressed archive](https://github.com/bioinformatics-centre/kaiju/archive/master.tar.gz) or
using the git command line client:
```
git clone https://github.com/bioinformatics-centre/kaiju.git
```
This will create the directory `kaiju` in the current directory.

Kaiju is written in C/C++11 for Linux and does not depend on additional libraries.
For compiling Kaiju and its associated programs, type:
```
cd kaiju/src
make
```
Afterwards, Kaiju's executable files are available in the `kaiju/bin` directory.
You can add this directory to your shell's `$PATH` variable or copy the files to a directory in your PATH.

##Creating the reference database and index

Before classification of reads, Kaiju's database index needs to be built from
the reference protein database.  You can either create a local index based on
the currently available data from GenBank, or download one of the indexes used
by the [Kaiju web server](http://kaiju.binf.ku.dk/).

For creating a local index, the program `makeDB.sh` in the `bin/` directory
will download the reference genomes and taxonomy files from the NCBI FTP server,
convert them into a protein database and construct Kaiju's index (the
Burrows-Wheeler transform and the FM-index) in one go.

It is recommended to create a new directory for downloading the files
and database construction, for example:
```
mkdir kaijudb
cd kaijudb
makeDB.sh [-r|-p|-n|-e]
```
The downloaded files are several GB in size. Therefore, the program should be
run in a directory having at least 80 GB of free space.

There are several options for creating the reference database with protein
sequences from different source databases:

###1. Complete Reference Genomes from NCBI RefSeq
`makeDB.sh -r`  
Download only completely assembled and annotated reference
genomes of Archaea and Bacteria from the NCBI RefSeq database.

Additionally, viral genomes from NCBI RefSeq can be added by using the option `-v`.

As of October 2016, this database contains ca. 20M protein sequences, which
amounts to a requirement of 14GB RAM for running Kaiju.

###2. Representative genomes from proGenomes
`makeDB.sh -p`  
Download the protein sequences belonging to the representative set of genomes
from the [proGenomes](http://progenomes.embl.de/) database.
This dataset generally covers a broader phylogenenic range compared to the RefSeq dataset,
and is therefore recommended, especially for environmental samples.

Additionally, viral genomes from NCBI RefSeq can be added by using the option `-v`.

As of October 2016, this database contains ca. 19M protein sequences, which
amounts to a requirement of 13GB RAM for running Kaiju.

###3. Non-redundant protein database _nr_
`makeDB.sh -n`  
Download the _nr_ database that is used by NCBI BLAST and extract proteins belonging
to Archaea, Bacteria and Viruses.

`makeDB.sh -e`  
Download the _nr_ database as above, but additionally include proteins from fungi and microbial eukaryotes.
The complete taxon list for this option is in the file `bin/taxonlist.tsv`.

Because the _nr_ database contains more proteins, more RAM is needed for index
construction and for running Kaiju.  As of October 2016, the _nr_ database with
option `-e` contains ca. 80M protein sequences, which amounts to a requirement
of 43GB RAM for running Kaiju.

###Index construction

When using option `-r`, `makeDB.sh` downloads and extracts 5 genomes from the NCBI FTP
server in parallel. This number can be changed by modifying the appropriate
variables at the beginning of the script.

By default, `makeDB.sh` uses 5 parallel threads for constructing the index, which can
be changed by using the option `-t`. Note that a higher number of threads
increases the memory usage during index construction, while reducing the number
of threads decreases memory usage.

After `makeDB.sh` is finished, only the files `kaiju_db.fmi` (or `kaiju_db_nr.fmi` / `kaiju_db_nr_euk.fmi`), `nodes.dmp`,
and `names.dmp` are needed to run Kaiju. The remaining files and the `genomes/`
directory containing the downloaded genomes can be deleted.

###Custom database
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
Then, Kaiju's index is created using the programs `mkbwt` and `mkfmi`. For example, if the database FASTA file is called `proteins.faa`, then run:
```
mkbwt -n 5 -a ACDEFGHIKLMNPQRSTVWY -o proteins proteins.faa
mkfmi proteins
```
which creates the file proteins.fmi that is used by Kaiju.
Note that the protein sequences may only contain the uppercase characters of the standard 20 amino acids, all other
characters need to be removed.

##Running Kaiju
Kaiju requires at least three arguments:
```
kaiju -t nodes.dmp -f kaiju_db.fmi -i inputfile.fastq
```
If you chose options `-n` or `-e` in `makeDB.sh`, then use `-f kaiju_db_nr.fmi` or `-f kaiju_db_nr_euk.fmi`.

For paired-end reads use `-i firstfile.fastq` and `-j secondfile.fastq`.

The reads must be in the same order in both files. Kaiju will strip suffixes
from the read names by deleting all characters after a `/` or space.  The read
names are then compared between the first and second file and an error is
issued if they are not identical.

Kaiju can read input files in FASTQ and FASTA format.
If the files are compressed, the shell's process substitution can be used to decompress them on the fly.
For example for GZIP compressed files:
```
kaiju -i <(gunzip -c firstfile.fastq.gz) -j <(gunzip -c secondfile.fastq.gz) ...
```

By default, Kaiju will print the output to the terminal (STDOUT).
The output can also be written to a file using the `-o` option:
```
kaiju -t nodes.dmp -f kaiju_db.fmi -i inputfile.fastq -o kaiju.out
```

Kaiju can use multiple parallel threads, which can be specified with the `-z` option, e.g. for using 25 parallel threads:
```
kaiju -z 25 -t nodes.dmp -f kaiju_db.fmi -i inputfile.fastq -o kaiju.out
```

###Run modes
The default run mode is **MEM**, which only considers exact matches.
For using the **Greedy** mode, which allows mismatches, set the mode via the option `-a` and the number
of allowed substitutions using option `-e`:
```
kaiju -t nodes.dmp -f kaiju_db.fmi -i inputfile.fastq -a greedy -e 5
```
The cutoffs for minimum required match length and match score can be changed using the options `-m` (default: 11) and `-s` (default: 65).

If the input sequences are already protein sequences, use option `-p` to disable translation of the input.

Option `-x` can be used to enable filtering of query sequences containing
low-complexity regions by using the SEG algorithm from the blast+ package.
Enabling this option is always recommended in order to avoid false positive
matches caused by spurious matches due to simple repeat patterns or other
sequencing noise.

###Output format
Kaiju will print one line for each read or read pair.
The default output format contains three columns separated by tabs.
Using the option `-v` enables the verbose output, which will print additional three columns:

1. either C or U, indicating whether the read is classified or unclassified.
2. name of the read
3. NCBI taxon identifier of the assigned taxon
4. the length or score of the best match used for classification
5. the taxon identifiers of all database sequences with the best match
6. matching fragment sequence(s)

##Classification accuracy

The accuracy of the classification depends both on the choice of the reference
database and the chosen options when running Kaiju. These choices also affect
the speed and memory usage of Kaiju.

For highest sensitivity, it is recommended to use the _nr_ database (+eukaryotes)
as a reference database because it is the most comprehensive set of protein
sequences. Alternatively, use proGenomes over Refseq for increased sensitivity.

Additionally, Greedy run mode, for example, with 5 allowed mismatches,
yields a higher sensitivity compared with MEM mode.

For fastest classification, use MEM mode and multiple parallel threads
(`-z`); and for lowest memory usage use the proGenomes reference
database. The number of parallel threads has only little impact on memory usage.

Further, the choice of the minimum required match length (`-m`) in MEM mode or
match score (`-s`) in Greedy mode governs the trade-off between sensitivity and
precision of the classification. Please refer to the paper for a discussion on
this topic.

##Helper programs
###Creating input file for Krona
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

###Creating classification summary
The program `kaijuReport` can convert Kaiju's tab-separated output file into a
summary report file for a given taxonomic rank, e.g., genus. It requires the
`nodes.dmp` and `names.dmp` files for mapping the taxon identifiers from
Kaiju's output to the corresponding taxon names.
```
kaijuReport -t nodes.dmp -n names.dmp -i kaiju.out -r genus -o kaiju.out.summary
```
The program can also filter out taxa with low abundances, e.g. for only showing genera that
comprise at least 1 percent of the total reads:
```
kaijuReport -t nodes.dmp -n names.dmp -i kaiju.out -r genus -m 1 -o kaiju.out.summary
```
or for showing genera comprising at least 1 percent of all classified reads:
```
kaijuReport -t nodes.dmp -n names.dmp -i kaiju.out -r genus -m 1 -u -o kaiju.out.summary
```
Option `-p` will print the full taxon path instead of just the taxon name.

###Adding taxa names
The program `addTaxonNames` appends the name that corresponds to the taxon id in
Kaiju's output file as a last column to the output.
```
addTaxonNames -t nodes.dmp -n names.dmp -i kaiju.out -o kaiju.names.out
```
Option `-u` will omit unclassified reads.  
Option `-p` will print the full taxon path instead of just the taxon name.  
Option `-r` will print the path containing only to the specified ranks. For example,
`-r phylum,genus` will append the names of phylum and genus to the end of each line.

###Merging outputs
The program `mergeOutputs` can merge two tab-separated output files in the
column format (see above) used by Kaiju and Kraken. Only the first three columns are used.

The files need to be sorted by the read name in the second column, for example by:
```
sort -k2,2 kaiju.out  >kaiju.out.sort
sort -k2,2 kraken.out >kraken.out.sort
```
Then both files can be merged:
```
mergeOutputs -i kaiju.out.sort -j kraken.out.sort -o combined.out -v
```
Again, process substitution can be used for sorting without creating intermediate files:
```
mergeOutputs -i <(sort -k2,2 kaiju.out) -j <(sort -k2,2 kraken.out) -o combined.out -v
```
The output file will be in the same column format as the input files (but only
contain the first three columns) and it will have the same length as the input
files (which also have to be of same length).  In the case of conflicting taxon
identifiers for a classified read in both input files, `mergeOutputs` will use the identifier found in the
first input file (specified by `-i`).  This behaviour can be changed using the
`-c` option, which can take four possible values:
- `1`: use taxon identifier from the first input file (default)
- `2`: use taxon identifier from the second input file
- `lca`: use the least common ancestor of the taxon identifiers from both files.
- `lowest`: use the lowest ranking of the two taxon identifiers iff they are within the same lineage. Otherwise use the LCA.

Options `lca` and `lowest` require the path to the file `nodes.dmp` by using the `-t` option.

###KaijuX and KaijuP

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
programs `mkbwt` and `mkfmi`:
```
mkbwt -n 5 -a ACDEFGHIKLMNPQRSTVWY -o proteins proteins.faa
mkfmi proteins
```
This will create two intermediate files `proteins.bwt` and `proteins.sa`, and finally
the file `proteins.fmi`, which is used by Kaiju.

The option `-n` for `mkbwt` specifies the number of parallel threads. The more
threads are used, the higher the memory consumption becomes.  The option `-e`
for `mkbwt` specifies the exponent of the suffix array checkpoint distances and
therefore determines the trade-off between the size of the suffix array and the
speed of the search. The default value is 5.

