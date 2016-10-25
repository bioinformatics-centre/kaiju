#!/bin/sh
SCRIPTDIR="$(dirname "$(readlink -f "$0")")"

viruses=0
use_nr=0
euk=0
use_progenomes=0
threadsBWT=5
parallelDL=5
parallelConversions=5
exponentSA=3
exponentSA_NR=5

usage() {
	echo This program creates a protein reference database and index for Kaiju.
	echo Several source databases can be used.
	echo
	echo By default, all complete bacterial and archaeal genomes in the
	echo NCBI RefSeq database are downloaded from the NCBI FTP server.
	echo
	echo Alternatively, all proteins belonging to the set of representative genomes
	echo from the proGenomes database can be downloaded using option -p.
  echo
  echo Viral proteins from NCBI Refseq can be included by using option -v.
	echo
	echo Instead of using full genomes as reference database, the NCBI BLAST
	echo non-redundant protein database \"nr\" can be downloaded using option -n.
	echo This file is then reduced to bacterial, archaeal, and viral proteins for
	echo Kaiju\'s database. Additionally, proteins belonging to fungi and microbial
	echo eukaryotes can also be included by using option -e.
	echo
	echo By default, 5 parallel threads are used for downloading and for the index construction.
	echo Use option -t to modify the number of threads. The more threads are used, the
	echo higher the memory requirement becomes.
	echo
	echo Optional arguments are:
	echo "  -p|--progenomes  download proteins from proGenomes database"
	echo "  -v|--viruses     download viral genomes from NCBI RefSeq"
	echo "  -n|--nr          download non-redundant protein database nr"
	echo "  -e|--euk         like -n, but also include microbial eukaryotes"
	echo "                   (listed in the file taxonlist.tsv)"
	echo "  -t|--threads X   for using X parallel threads for index construction"
	echo
}

while :; do
    case $1 in
        -h|-\?|--help)
            usage
            exit 1
            ;;
        -t|--threads)
            if [ -n "$2" ]; then
                threadsBWT=$2
                shift
            else
                printf 'ERROR: "--threads" requires a non-empty integer argument.\n' >&2
                usage
                exit 1
            fi
            ;;
        -n|--nr)
            use_nr=1
            ;;
        -e|--euk)
            euk=1
            ;;
        -v|--viruses)
            viruses=1
            ;;
        -p|--progenomes)
            use_progenomes=1
            ;;
        --)# End of all options.
            shift
            break
            ;;
        -?*)
            printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
            ;;
        *)# Default case: If no more options then break out of the loop.
            break
    esac
    shift
done

#check if necessary programs are in the PATH
command -v awk >/dev/null 2>/dev/null || { echo Error: awk not found; exit 1; }
command -v wget >/dev/null 2>/dev/null || { echo Error: wget not found; exit 1; }
command -v xargs >/dev/null 2>/dev/null || { echo Error: xargs not found; exit 1; }
command -v tar >/dev/null 2>/dev/null || { echo Error: tar not found; exit 1; }
command -v gunzip >/dev/null 2>/dev/null || { echo Error: gunzip not found; exit 1; }
command -v perl >/dev/null 2>/dev/null || { echo Error: perl not found; exit 1; }

command -v $SCRIPTDIR/gbk2faa.pl >/dev/null 2>/dev/null || { echo Error: gbk2faa.pl not found in $SCRIPTDIR; exit 1; }
command -v $SCRIPTDIR/mkfmi >/dev/null 2>/dev/null || { echo Error: mkfmi not found in $SCRIPTDIR; exit 1; }
command -v $SCRIPTDIR/mkbwt >/dev/null 2>/dev/null || { echo Error: mkbwt not found in $SCRIPTDIR; exit 1; }
command -v $SCRIPTDIR/convertNR >/dev/null 2>/dev/null || { echo Error: convertNR not found in $SCRIPTDIR; exit 1; }

if [ "$euk" -eq 1 ]
then
	[ -r $SCRIPTDIR/taxonlist.tsv ] || { echo Error: file taxonlist.tsv not found in $SCRIPTDIR; exit 1; }
fi

#test AnyUncompress usable in perl
`perl -e 'use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);'`
[ $? -ne 0 ] && { echo Error: Perl IO::Uncompress::AnyUncompress library not found; exit 1; }

#good to go
set -e

echo Downloading file: taxdump.tar.gz
wget --show-progress -N -nv ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
if [ -r taxdump.tar.gz ]
then
	echo Extracting nodes.dmp and names.dmp files
	tar xf taxdump.tar.gz nodes.dmp names.dmp
else 
	echo Error: File taxdump.tar.gz was not downloaded
	exit 1
fi

if [ "$use_nr" -eq 1 -o "$euk" -eq 1 ]
then
	echo Downloading file nr.gz
	wget --show-progress -N -c -nv ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
	echo Downloading file prot.accession2taxid.gz
	wget --show-progress -N -c -nv ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
	if [ -r nr.gz -a -r prot.accession2taxid.gz ]
	then
		echo Unpacking prot.accession2taxid.gz
		gunzip -f prot.accession2taxid.gz
		echo Converting NR file to Kaiju database
		if [ "$euk" -eq 1 ]
		then
			gunzip -c nr.gz | $SCRIPTDIR/convertNR -t nodes.dmp -g prot.accession2taxid -c -o kaiju_db_nr_euk.faa -l $SCRIPTDIR/taxonlist.tsv
			echo Creating BWT from Kaiju database
			$SCRIPTDIR/mkbwt -e $exponentSA_NR -n $threadsBWT -a ACDEFGHIKLMNPQRSTVWY -o kaiju_db_nr_euk kaiju_db_nr_euk.faa
			echo Creating FM-index
			$SCRIPTDIR/mkfmi kaiju_db_nr_euk
			echo Done!
			echo You can delete the files nr.gz, taxdump.tar.gz, prot.accession2taxid, kaiju_db_nr_euk.faa, kaiju_db_nr_euk.bwt, kaiju_db_nr_euk.sa
			echo Kaiju only needs the files kaiju_db_nr_euk.fmi, nodes.dmp, and names.dmp.
		else
			gunzip -c nr.gz | $SCRIPTDIR/convertNR -t nodes.dmp -g prot.accession2taxid -c -o kaiju_db_nr.faa
			echo Creating BWT from Kaiju database
			$SCRIPTDIR/mkbwt -e $exponentSA_NR -n $threadsBWT -a ACDEFGHIKLMNPQRSTVWY -o kaiju_db_nr kaiju_db_nr.faa
			echo Creating FM-index
			$SCRIPTDIR/mkfmi kaiju_db_nr
			echo Done!
			echo You can delete the files nr.gz, taxdump.tar.gz, prot.accession2taxid, kaiju_db_nr.faa, kaiju_db_nr.bwt, kaiju_db_nr.sa
			echo Kaiju only needs the files kaiju_db_nr.fmi, nodes.dmp, and names.dmp.
		fi
	fi
else
	echo Creating directory genomes/
	mkdir -p genomes
	if [ "$use_progenomes" -eq 0 ]
	then
		echo Downloading file list for complete genomes from RefSeq...
		wget -nv -O assembly_summary.archaea.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt
		wget -nv -O assembly_summary.bacteria.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
		awk 'BEGIN{FS="\t";OFS="/"}$12=="Complete Genome" && $11=="latest"{l=split($20,a,"/");print $20,a[l]"_genomic.gbff.gz"}' assembly_summary.bacteria.txt assembly_summary.archaea.txt > downloadlist.txt
		nfiles=`cat downloadlist.txt| wc -l`
		echo Downloading $nfiles genome files from GenBank FTP server. This may take a while...
		cat downloadlist.txt | xargs -P $parallelDL -n 1 wget -P genomes -c -nv

		if [ "$viruses" -eq 1 ]
		then
			echo Downloading virus genomes from RefSeq...
			wget --show-progress -N -c -nv -P genomes ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.genomic.gbff.gz
			wget --show-progress -N -c -nv -P genomes ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.genomic.gbff.gz
		fi

		echo Extracting protein sequences from downloaded files...
		find ./genomes -name "*.gbff.gz" | xargs -n 1 -P $parallelConversions -i $SCRIPTDIR/gbk2faa.pl '{}' '{}'.faa
	else
		echo Downloading proGenomes database...
		wget --show-progress -N -c -nv -P genomes http://progenomes.embl.de/data/repGenomes/representatives.proteins.fasta.gz

		if [ "$viruses" -eq 1 ]
		then
			echo Downloading virus genomes from RefSeq...
			wget --show-progress -N -c -nv -P genomes ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.genomic.gbff.gz
			wget --show-progress -N -c -nv -P genomes ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.genomic.gbff.gz
		fi

		echo Extracting protein sequences from downloaded files...
		gunzip -c genomes/representatives.proteins.fasta.gz | perl -lsane 'if(/>(\d+)\./){print ">",++$c,"_",$1}else{y/BZ/DE/;s/[^ARNDCQEGHILKMFPSTWYV]//gi;print if length}' > genomes/representatives.proteins.fasta.gz.faa
		find ./genomes -name "viral.*.gbff.gz" | xargs -n 1 -P $parallelConversions -i $SCRIPTDIR/gbk2faa.pl '{}' '{}'.faa
	fi

	cat genomes/*.faa >kaiju_db.faa

	echo Creating Borrows-Wheeler transform...
	$SCRIPTDIR/mkbwt -n $threadsBWT -e $exponentSA -a ACDEFGHIKLMNPQRSTVWY -o kaiju_db kaiju_db.faa
	echo Creating FM-Index...
	$SCRIPTDIR/mkfmi kaiju_db
	echo Done!
	echo You can delete the folder genomes/ as well as the files taxdump.tar.gz, kaiju_db.faa, kaiju_db.bwt, and kaiju_db.sa
	echo Kaiju only needs the files kaiju_db.fmi, nodes.dmp, and names.dmp.
fi


