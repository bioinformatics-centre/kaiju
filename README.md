#Kaiju

Authors:  
Peter Menzel <pmenzel@gmail.com>   
Anders Krogh <krogh@binf.ku.dk>   


Kaiju is a program for the taxonomic assignment of high-throughput sequencing
reads, e.g., Illumina or Roche/454, from whole-genome sequencing of
metagenomic DNA. Reads are directly assigned to taxa using the NCBI taxonomy and a 
reference database of protein sequences from bacterial and archaeal genomes.

The program is described in a [preprint paper available on bioRxiv](http://biorxiv.org/content/early/2015/11/16/031229).

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

Kaiju is written in C/C++ for Linux and does not depend on additional libraries.
For compiling Kaiju and its associated programs, type:
```
cd kaiju/src
make
```
Afterwards, Kaiju's executable files will located in the `kaiju/bin` directory.
You can add this directory to your shell's `$PATH` variable or copy the files to a directory in your PATH.

##Creating the reference database and index

Before classification of reads, Kaiju's database index needs to be built from
the reference protein database.  You can either create a local index based on
the currently available data from GenBank, or download the current index used
by the [Kaiju web server](http://kaiju.binf.ku.dk/).

For creating a local index, the program `makeDB.sh` in the `bin/` directory
will download the complete genome and taxonomy files from the NCBI FTP server,
convert them to the protein database and construct Kaiju's index (the
Burrows-Wheeler-transform and the FM-index) in one go.

It is recommended to create a new folder for the download and run the program from there, e.g.:
```
mkdir kaijudb
cd kaijudb
makeDB.sh
```
The downloaded files are several GB in size. Therefore, the program should be
run in a directory with at least 50 GB free space.

There are two options for the reference database:
###1. Complete Genomes
The first option is to use only completely assembled and annotated reference genomes
from the RefSeq database in GenBank. This is the default behaviour of `makeDB.sh`.
Additional to archaeal and bacterial genomes, `makeDB.sh` can also add viral genomes to the database by using the option `-v`.
As of Feb 2015, this database contains ca. 15m protein sequences, which amounts to ca. 9GB memory
required by Kaiju.

By default, `makeDB` downloads and extracts 5 files in parallel. This number can
be changed by modifying the appropriate variables at the beginning of the
script.  The program also uses 5 parallel threads for construction the index,
which can be changed by using the option `-t`. Note that a higher number of threads
increases the memory usage during index construction.

After `makeDB.sh` is finished, only the files `kaiju_db.fmi`, `nodes.dmp`,
and `names.dmp` are needed to run Kaiju.  The remaining files and the `genomes`
folder containing the downloaded genomes can be deleted.

###2. Non-redundant protein database
The second option is to use the complete protein database (NR) from GenBank
as a refererence database. Use the option `-n` for `makeDB.sh`.
The program will download the `nr.gz` file from GenBank's FTP server and convert it into
a database by excluding all proteins that are not assigned to Bacteria, Archaea, or Viruses
in the NCBI taxonomy. 
Since NR contains many more proteins, more memory is
needed for index construction and for running Kaiju.
As of Feb 2015, this database contains ca. 61m protein sequences, which amounts to ca. 31GB memory
required by Kaiju.

After `makeDB.sh` is finished, only the files `kaiju_db_nr.fmi`, `nodes.dmp`,
and `names.dmp` are needed to run Kaiju.  The remaining files and the `genomes`
folder containing the downloaded genomes can be deleted.


##Running Kaiju
Kaiju requires at least three arguments:
```
kaiju -t nodes.dmp -f kaiju_db.fmi -i inputfile.fastq
```
If you choose to use the NR database, then use `-f kaiju_db_nr.fmi`.

For paired-end reads use `-i firstfile.fastq` and `-j secondfile.fastq`.
The reads must be in the same order in both files.

Kaiju can read input files in FASTQ and FASTA format.
If the files are compressed, the shell's process substitution can be used to decompress them on the fly.
For example for GZIP compressed files:
```
kaiju ... -i <(gunzip -c firstfile.fastq.gz) -j <(gunzip -c secondfile.fastq.gz) ...
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

The default run mode is MEM. For using Greedy mode, set the mode via the option `-a` and the number
of allowed substitutions using option `-e`:
```
kaiju -t nodes.dmp -f kaiju_db.fmi -i inputfile.fastq -a greedy -e 5
```
The cutoffs for minimum required match length and match score can be changed using the options `-m` and `-s`.



###Output format
Kaiju will print one line for each read or read pair.
The default output format contains three columns separated by tabs.
Using the option `-v` enables the verbose output, which will print additional three columns:

1. either C or U, indicating whether the read is classified or unclassified.
2. name of the read
3. NCBI taxon identifier of the assigned taxon
4. the length or score of the best match used for classification
5. the taxon identifiers of all database sequences with the best match
6. matching fragment(s) sequence

##Helper programs
###Creating input file for Krona
The program `kaiju2krona` can be used to convert Kaiju's tab-separated output file
into a tab-separated text file, which can be imported into [Krona](https://github.com/marbl/Krona/wiki/KronaTools). It requires the nodes.dmp
and names.dmp files from the NCBI taxonomy in order to map the taxon identifiers from Kaiju's
output to the taxon names.
```
kaiju2krona -t nodes.dmp -n names.dmp -i kaiju.out -o kaiju.out.krona
```
The file `kaiju.out.krona` can then be imported into Krona and converted into an HTML file using
Krona's `ktImportText` program:
```
ktImportText -o kaiju.out.html kaiju.out.krona
```

###Creating summary
The program `kaijuReport` can convert Kaiju's tab-separated output file
into a summary report file for a given taxonomic rank, e.g., genus. It requires the nodes.dmp
and names.dmp files from the NCBI taxonomy in order to map the taxon identifiers from Kaiju's
output to the taxon names.
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
Again, process substitution can be used without creating intermediate files:
```
mergeOutputs -i <(sort -k2,2 kaiju.out) -j <(sort -k2,2 kraken.out) -o combined.out -v
```
The output file will be in the same column format as the input files (but only
contain the first three columns) and it will have the same length as the input
files (which have to be of same length).  In the case of conflicting taxon identifiers in both files,
`mergeOutputs` will by default use the identifier found in the first input file (specified by `-i`).
This behavior can be changed by the `-c` option, which can take the values
`1` (default), `2` (use identifier from the second file) or `lca`, which determines and prints
the least common ancestor of the taxon identifiers from both files. Option `lca`
requires to specify the nodes.dmp file using the `-t` option.

###KaijuX and KaijuP

The programs `kaijux` and `kaijup` can be used for finding the best matches for a query sequence
in a protein database without taxonomic classification, i.e., they will just print the identifier
of the database sequence. Thus, these programs do not use the nodes.dmp file containing the taxonomy,
but only need the `.fmi` database file. While `kaijux` takes nucleotide sequences as input and translates
them into the six reading frames, just like standard `kaiju`, `kaijup` takes protein sequences as input,
which are directly searched in the database.
All other parameters remain the same as in standard `kaiju`. In case of paired-end reads, both mates are
searched independently.


