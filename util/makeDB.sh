#!/bin/sh
SCRIPTDIR="$(dirname "$(readlink -f "$0")")"

viruses=0
threadsBWT=5
parallelDL=5
parallelConversions=5
exponentSA=3

usage() {
	echo This programs downloads the complete bacterial and archaeal genomes from
	echo the GenBank FTP server and makes a database index for Kaiju.
  echo By default, it uses 5 parallel threads for downloading and for index construction.
  echo Optional arguments are:
  echo   -t\|--threads X   for using X parallel threads for index construction.
  echo   -v\|--viruses     for also downloading viral genomes  
}

while :; do
    case $1 in
        -h|-\?|--help)   # Call a "show_help" function to display a synopsis, then exit.
            usage
            exit
            ;;
        -t|--threads)       # Takes an option argument, ensuring it has been specified.
            if [ -n "$2" ]; then
                threadsBWT=$2
                shift
            else
                printf 'ERROR: "--threads" requires a non-empty integer argument.\n' >&2
								usage
                exit 1
            fi
            ;;
        -v|--viruses)
            viruses=1
            ;;
        --)              # End of all options.
            shift
            break
            ;;
        -?*)
            printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
            ;;
        *)               # Default case: If no more options then break out of the loop.
            break
    esac
    shift
done

#check if necessary programs are in the PATH
command -v awk >/dev/null 2>/dev/null || { echo Error: awk not found; exit;}
command -v wget >/dev/null 2>/dev/null || { echo Error: wget not found; exit;}
command -v xargs >/dev/null 2>/dev/null || { echo Error: xargs not found; exit;}

command -v $SCRIPTDIR/gbk2faa.pl >/dev/null 2>/dev/null || { echo Error: gbk2faa.pl not found in $SCRIPTDIR; exit;}
command -v $SCRIPTDIR/mkfmi >/dev/null 2>/dev/null || { echo Error: mkfmi not found in $SCRIPTDIR; exit;}
command -v $SCRIPTDIR/mkbwt >/dev/null 2>/dev/null || { echo Error: mkbwt not found in $SCRIPTDIR; exit;}

#test AnyUncompress usable in perl
`perl -e 'use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);'`
if [ $? -ne 0 ]
then
	echo Error: Perl IO::Uncompress::AnyUncompress library not found
	exit
fi

#good to go
set -e

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
$SCRIPTDIR/mkbwt -n $threadsBWT -e $exponentSA -a ACDEFGHIKLMNPQRSTVWY -o allproteins allproteins.faa
echo Creating FM-Index...
$SCRIPTDIR/mkfmi  allproteins

echo Downloading taxonomy files from NCBI
wget -nv ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
if [ -r taxdump.tar.gz ]
then
	echo Extracting nodes.dmp and names.dmp files
	tar xf taxdump.tar.gz nodes.dmp names.dmp
else 
	echo Error: File taxdump.tar.gz was not downloaded
	exit
fi

echo Done!
echo You can delete the genomes folder as well as allproteins.faa, allproteins.bwt, and allproteins.sa
echo Kaiju only needs the files allproteins.fmi, nodes.dmp, and names.dmp.

