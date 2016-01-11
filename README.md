
#Kaiju

Authors:  
Peter Menzel <pmenzel@gmail.com>   
Anders Krogh <krogh@binf.ku.dk>   


Kaiju is a program for the taxonomic assignment of high-throughput sequencing
reads, e.g., Illumina or Roche/454, from whole-genome sequencing of
metagenomic DNA. Reads are directly assigned to taxa using the NCBI taxonomy and a 
reference database of protein sequences from bacterial and archaeal genomes.

The program is described in a [preprint paper available on bioRxiv](http://biorxiv.org/content/early/2015/11/16/031229).

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
Kaiju can be downloaded directly from GitHub either as a compressed archive or 
using the git command line client:
```
git clone https://github.com/bioinformatics-centre/kaiju.git
```
This will create the directory `kaiju` in the current directory.

Kaiju is written in C and C++ for Linux and does not depend on additional libraries.   
For compiling Kaiju and its associated programs, type:
```
cd kaiju/src
make
```
This will create the executable files `mkfmi`, `mkbwt`, `kaiju`, `kaijux`, `kaijup`, `mergeOutputs`, `kaijuReport`, 
`kaiju2krona`, and `gbk2faa.pl` in the `kaiju/bin` directory.
You can add this directory to your shell's PATH variable or copy the files to a directory in your PATH.

##Creating the reference database
###Download files
The first step is to download the GenBank files for microbial genomes
and the taxonomy information from the NCBI FTP server. 
The files are several GB in size and need to be decompressed. Therefore, these steps
need to be done in a directory with at least 50 GB free space.
```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Bacteria/all.gbk.tar.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
```
Note that the location of the files on the server might change in the future.

Then the GenBank files are extracted into a subdirectory `genomes` and 
the needed taxonomy files are extracted to the current directory:
```
mkdir genomes
tar -C genomes -xf all.gbk.tar.gz
tar xf taxdump.tar.gz nodes.dmp names.dmp
```
###Convert gbk files
The Perl script `gbk2faa.pl` converts
`.gbk` files into `.faa` files, which have the taxon id in the FASTA header lines.

Run this script for each of the GenBank files and finally paste all sequences
into one file `allproteins.faa`:

```
find ./genomes -name "*.gbk" | xargs -i gbk2faa.pl '{}' '{}'.faa
cat all/*/*.faa >allproteins.faa
```
Several files can be converted simultaneously by using the option `-P` for `xargs`.

###Create Borrows-Wheeler-transform and FM-index
The last step creates the database index (FM-index and suffix array) that is used by Kaiju from `allproteins.faa`:
```
mkbwt -e 3 -n 5 -a ACDEFGHIKLMNPQRSTVWY -o allproteins allproteins.faa
mkfmi -i allproteins
```
The `-n` option of `mkbwt` specifies the number of parallel CPU threads to use. The more threads are used, the higher the memory consumption is.
For example, 5 parallel threads require around 7.5 GB memory with the current NCBI database comprising 2785 genomes.

The `-e` option specifies the exponent of the suffix array checkpoint distances and therefore influences the size of the suffix array and
the speed of the search. A value of 3 translates roughly to a database size of ca. 5 GB and increasing the number will decrease the size.

The file allproteins.fmi created by `mkfmi` and nodes.dmp are the only files needed to run Kaiju.
The remaining protein files can be deleted.

##Running Kaiju
Kaiju requires at least three arguments:
```
kaiju -t nodes.dmp -f allproteins.fmi -i inputfile.fastq
```
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
kaiju -t nodes.dmp -f allproteins.fmi -i inputfile.fastq -o kaiju.out
```

The default run mode is MEM. For using Greedy mode, set the mode via the option `-a` and the number 
of allowed substitutions using option `-e`:
```
kaiju -t nodes.dmp -f allproteins.fmi -i inputfile.fastq -a greedy -e 5
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


