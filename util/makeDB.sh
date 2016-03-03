#!/bin/sh
SCRIPTDIR="$(dirname "$(readlink -f "$0")")"

viruses=0
use_nr=0
threadsBWT=5
parallelDL=5
parallelConversions=5
exponentSA=3
exponentSA_NR=5

usage() {
	echo This program downloads all complete bacterial and archaeal genomes from
	echo the GenBank FTP server and builds a database index for Kaiju.
	echo
	echo Alternatively, makeDB.sh can download the non-redundant protein database NR
	echo from GenBank, which contains all proteins. This file is then reduced to
	echo bacterial, archaeal, and viral proteins for Kaiju\'s database.
	echo
	echo By default, this program uses 5 parallel threads for downloading and for
	echo index construction.
	echo
	echo Optional arguments are:
	echo "  -t|--threads X   for using X parallel threads for index construction"
	echo "  -v|--viruses     for also downloading viral genomes"
	echo "  -n|--nr          download GenBank non-redundant protein database instead"
}

while :; do
    case $1 in
        -h|-\?|--help)
            usage
            exit
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
        -v|--viruses)
            viruses=1
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
command -v awk >/dev/null 2>/dev/null || { echo Error: awk not found; exit;}
command -v wget >/dev/null 2>/dev/null || { echo Error: wget not found; exit;}
command -v xargs >/dev/null 2>/dev/null || { echo Error: xargs not found; exit;}
command -v tar >/dev/null 2>/dev/null || { echo Error: tar not found; exit;}
command -v gunzip >/dev/null 2>/dev/null || { echo Error: gunzip not found; exit;}
command -v perl >/dev/null 2>/dev/null || { echo Error: perl not found; exit;}

command -v $SCRIPTDIR/gbk2faa.pl >/dev/null 2>/dev/null || { echo Error: gbk2faa.pl not found in $SCRIPTDIR; exit;}
command -v $SCRIPTDIR/mkfmi >/dev/null 2>/dev/null || { echo Error: mkfmi not found in $SCRIPTDIR; exit;}
command -v $SCRIPTDIR/mkbwt >/dev/null 2>/dev/null || { echo Error: mkbwt not found in $SCRIPTDIR; exit;}
command -v $SCRIPTDIR/convertNR >/dev/null 2>/dev/null || { echo Error: convertNR not found in $SCRIPTDIR; exit;}

#test AnyUncompress usable in perl
`perl -e 'use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);'`
if [ $? -ne 0 ]
then
	echo Error: Perl IO::Uncompress::AnyUncompress library not found
	exit
fi

#good to go
set -e

echo Downloading taxonomy files from NCBI
wget -N -nv ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
if [ -r taxdump.tar.gz ]
then
	echo Extracting nodes.dmp and names.dmp files
	tar xf taxdump.tar.gz nodes.dmp names.dmp
else 
	echo Error: File taxdump.tar.gz was not downloaded
	exit
fi

if [ "$use_nr" -eq 1 ]
then
	echo Downloading NR file...
	wget -N -c -nv ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
	echo Downloading gi_taxid_prot.dmp file...
	wget -N -c -nv ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_prot.dmp.gz
	if [ -r nr.gz -a -r gi_taxid_prot.dmp.gz ]
	then
		echo Unpacking gi_taxid_prot.dmp.gz 
		gunzip gi_taxid_prot.dmp.gz
		echo Converting NR file to Kaiju database
		gunzip -c nr.gz | $SCRIPTDIR/convertNR -t nodes.dmp -g gi_taxid_prot.dmp -c -o kaiju_db_nr.faa
		echo Creating BWT from Kaiju database
		$SCRIPTDIR/mkbwt -e $exponentSA_NR -n $threadsBWT -a ACDEFGHIKLMNPQRSTVWY -o kaiju_db_nr kaiju_db_nr.faa
		echo Creating FM-index
		$SCRIPTDIR/mkfmi kaiju_db_nr
		echo Done!
		echo You can delete the files nr.gz, taxdump.tar.gz, gi_taxid_prot.dmp, kaiju_db_nr.faa, kaiju_db_nr.bwt, kaiju_db_nr.sa
		echo Kaiju only needs the files kaiju_db_nr.fmi, nodes.dmp, and names.dmp.
	fi
else
	echo Creating directory genomes/
	mkdir -p genomes
	echo Downloading file list for full genomes...
	wget -nv -O assembly_summary.archaea.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt
	wget -nv -O assembly_summary.bacteria.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
	awk 'BEGIN{FS="\t";OFS="/"}$12=="Complete Genome" && $11=="latest"{split($20,a,"/");print $20,a[6]"_genomic.gbff.gz"}' assembly_summary.bacteria.txt assembly_summary.archaea.txt > downloadlist.txt
	nfiles=`cat downloadlist.txt| wc -l`
	echo Downloading $nfiles genome files from GenBank FTP server. This may take a while...
	cat downloadlist.txt | xargs -P $parallelDL -n 1 wget -P genomes -c -nv

	if [ "$viruses" -eq 1 ]
	then
		echo Downloading virus genomes from GenBank FTP server...
		wget -c -nv -P genomes ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.genomic.gbff.gz
		wget -c -nv -P genomes ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.genomic.gbff.gz
	fi

	echo Extracting protein sequences from downloaded files...
	find ./genomes -name "*.gbff.gz" | xargs -n 1 -P $parallelConversions -i $SCRIPTDIR/gbk2faa.pl '{}' '{}'.faa
	cat genomes/*.gz.faa >allproteins.faa

	echo Creating Borrows-Wheeler transform...
	$SCRIPTDIR/mkbwt -n $threadsBWT -e $exponentSA -a ACDEFGHIKLMNPQRSTVWY -o kaiju_db allproteins.faa
	echo Creating FM-Index...
	$SCRIPTDIR/mkfmi kaiju_db
	echo Done!
	echo You can delete the folder genomes/ as well as the files allproteins.faa, kaiju_db.bwt, and kaiju_db.sa
	echo Kaiju only needs the files kaiju_db.fmi, nodes.dmp, and names.dmp.
fi


