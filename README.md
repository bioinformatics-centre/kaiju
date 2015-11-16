
#Kaiju

version 1.0

Authors:  
Peter Menzel <pmenzel@gmail.com>   
Anders Krogh <krogh@binf.ku.dk>   


Kaiju is a program for the taxonomic assignment of high-throughput sequencing
reads, e.g., Illumina or Roche/454, from whole-genome sequencing of
metagenomic DNA. Reads are directly assigned to taxa using the NCBI taxonomy and a 
reference database of protein sequences from bacterial and archaeal genomes.

###License

Copyright (c) 2015 Peter Menzel and Anders Krogh

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
For compiling Kaiju and its associated programs:
```
cd kaiju/src
make
```
This will create the executable files `mkfmi`, `mkbwt`, `kaiju`, `mergeOutputs`, and
`kaiju2krona` in the `kaiju/bin` directory.
You can add this directory to your shell's `PATH` variable or copy the those files to a directory in your PATH.

##Creating the reference database
###Download files
The first step is to download the GenBank files for microbial genomes
and the taxonomy information from the NCBI FTP server. 
The files are several GB in size and need to be decompressed. Therefore these steps
need to be done in a directory with at least 50GB free space.
```
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Bacteria/all.gbk.tar.gz
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
The Perl script `gbk2faa.pl` in the `kaiju/utils/` subdirectory converts
gbk files into faa files, which have the taxon id in the Fasta headers.

Run this script for each of the GenBank files and finally paste all sequences
into one file `allproteins.faa`:

```
find ./genomes -name "*.gbk" | xargs -i gbk2faa.pl '{}' '{}'.faa
cat all/*/*.faa >allproteins.faa
```

If GNU parallel is installed, it can be used for simultaneously converting several files, e.g., for 5 files at the same time:
```
find ./genomes -name "*.gbk" | parallel -j 5 gbk2faa.pl '{}' '{}'.faa
cat all/*/*.faa >allproteins.faa
```

###Create Borrows-Wheeler-transform and FM-index
The last step creates the index files (FM-index and suffix array) that are used by Kaiju from `allproteins.faa`:
```
mkbwt -n 20 -a ACDEFGHIKLMNPQRSTVWY allproteins.faa > allproteins.bwt
mkfmi -n 20 -i allproteins.bwt allproteins
```
The `-n` option specifies the number of parallel  CPU threads to use, for example 20.
The two files `allproteins.fmi` and `allproteins.sa` are created by  `mkfmi` and only those two are used by Kaiju.
The remaining protein files can be deleted.

##Running Kaiju
Kaiju can be run in default mode using these four mandatory options:
```
kaiju -t nodes.dmp -f allproteins.fmi -b allproteins.sa -i inputfile.fastq
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
kaiju -t nodes.dmp -f allproteins.fmi -b allproteins.sa -i inputfile.fastq -o kaiju.out
```

The default run mode is MEM. For using Greedy modes with mismatches set the mode via the option `-a` and the number 
of allowed substitutions using option `-e`:
```
kaiju -t nodes.dmp -f allproteins.fmi -b allproteins.sa -i inputfile.fastq -a greedy -e 5
```
The cutoffs for minimum required match length and match score can be changed using the options `-m` and `-s`.


###Output format
Kaiju will print one line for each read or read pair.
The default output format contains three columns.
Using the option `-v` enables the verbose output, which will print additional three columns:

1. either C or U, indicating whether the read is classified or unclassified. 
2. name of the read
3. NCBI taxon identifier of the assigned taxon
4. the length or score of the best match used for classification
5. the taxon identifiers of all database sequences with the best match
6. matching fragment(s) sequence

##Helper programs
###Creating input file for Krona
The program `kaiju2krona` can be used to convert Kaiju's output file
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

###Merging outputs
The program `mergeOutputs` can merge two output files of the same length in the
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
The output file will also be in the three column format and have the same
length as the input files.  In the case of conflicting taxon identifiers in both files,
`mergeOutputs` will by default use the identifier found in the first input file (specified by `-i`).
This behavior can be changed by the `-c` option, which can take the values
`1` (default), `2` (use identifier from the second file) or `lca`, which determines and prints
the least common ancestor of the taxon identifiers from both files. Option `lca`
requires to specify the nodes.dmp file using the `-t` option.

