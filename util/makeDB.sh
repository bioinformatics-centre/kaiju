#!/bin/sh
#
# This file is part of Kaiju, Copyright 2015-2017 Peter Menzel and Anders Krogh
# Kaiju is licensed under the GPLv3, see the file LICENSE.
#
SCRIPTDIR=$(dirname $0)

PATH=$PATH:$SCRIPTDIR

db_viruses=0
db_refseq=0
db_progenomes=0
db_nr=0
db_euk=0
threadsBWT=5
parallelDL=5
parallelConversions=5
exponentSA=3
exponentSA_NR=5
DL=1

usage() {
s=" "
tab="    "
echo
echo This program creates a protein reference database and index for Kaiju.
echo Several source databases can be used and one of these options must be set:
echo
echo  "$s" -r  all complete bacterial and archaeal genomes in the NCBI RefSeq database
echo
echo  "$s" -p  all proteins belonging to the set of representative genomes
echo  "$tab"   from the proGenomes database
echo
echo  "$s" -n  NCBI BLAST non-redundant protein database \"nr\":
echo  "$tab"   only Archaea, bacteria, and viruses
echo
echo  "$s" -e   NCBI BLAST non-redundant protein database \"nr\":
echo  "$tab"   like -n, but additionally including fungi and microbial eukaryotes
echo
echo Additional options:
echo
echo  "$s" -v    additionally add viral genomes from RefSeq,
echo  "$tab"   when using the RefSeq or proGenomes database
echo
echo  "$s" -t X  set number of parallel threads for index construction to X \(default:5\)
echo  "$tab"     The more threads are used, the higher the memory requirement becomes.
echo
echo  "$s" --noDL  do not download files, but use the existing files in the folder.
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
                printf 'ERROR: "-t" requires a non-empty integer argument.\n' >&2
                usage
                exit 1
            fi
            ;;
        --noDL)
            DL=0
            ;;
        -n|--nr)
            db_nr=1
            ;;
        -e|--euk)
            db_euk=1
            ;;
        -v|--viruses)
            db_viruses=1
            ;;
        -p|--progenomes)
            db_progenomes=1
            ;;
        -r|--refseq)
            db_refseq=1
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

[ $db_refseq -eq 1 -o $db_progenomes -eq 1 -o $db_nr -eq 1 -o $db_euk -eq 1 ] || { echo "Error: Use one of the options -r, -p, -n or -e"; usage; exit 1; }

#check if necessary programs are in the PATH
command -v awk >/dev/null 2>/dev/null || { echo Error: awk not found; exit 1; }
command -v wget >/dev/null 2>/dev/null || { echo Error: wget not found; exit 1; }
command -v xargs >/dev/null 2>/dev/null || { echo Error: xargs not found; exit 1; }
command -v tar >/dev/null 2>/dev/null || { echo Error: tar not found; exit 1; }
command -v gunzip >/dev/null 2>/dev/null || { echo Error: gunzip not found; exit 1; }
command -v perl >/dev/null 2>/dev/null || { echo Error: perl not found; exit 1; }

#test if option --show-progress is available for wget, then use it when downloading
wgetProgress=""
wget --help | grep -q -- --show-progress && wgetProgress='--show-progress'

#check that all programs from Kaiju are usable
command -v gbk2faa.pl >/dev/null 2>/dev/null || { echo Error: gbk2faa.pl not found in $PATH; exit 1; }
command -v mkfmi >/dev/null 2>/dev/null || { echo Error: mkfmi not found in $PATH; exit 1; }
command -v mkbwt >/dev/null 2>/dev/null || { echo Error: mkbwt not found in $PATH; exit 1; }
command -v convertNR >/dev/null 2>/dev/null || { echo Error: convertNR not found in $PATH; exit 1; }

if [ $db_euk -eq 1 ]
then
	[ -r $SCRIPTDIR/taxonlist.tsv ] || { echo Error: file taxonlist.tsv not found in $SCRIPTDIR; exit 1; }
fi

#test AnyUncompress usable in perl
`perl -e 'use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError);'`
[ $? -ne 0 ] && { echo Error: Perl IO::Uncompress::AnyUncompress library not found; exit 1; }

#good to go
set -e

if [ $DL -eq 1 ]
then
	echo Downloading file taxdump.tar.gz
	wget -N -nv $wgetProgress ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
fi
[ -r taxdump.tar.gz ] || { echo Missing file taxdump.tgz; exit 1; }
echo Extracting file taxdump.tar.gz
tar xf taxdump.tar.gz nodes.dmp names.dmp merged.dmp

if [ $db_nr -eq 1 -o $db_euk -eq 1 ]
then
	if [ $DL -eq 1 ]
	then
		echo Downloading file nr.gz
		wget -N -nv $wgetProgress ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
		echo Downloading file prot.accession2taxid.gz
		wget -N -nv $wgetProgress ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
	fi
	[ -r nr.gz ] || { echo Missing file nr.gz; exit 1; }
	[ -r prot.accession2taxid.gz ] || { echo Missing file prot.accession2taxid.gz; exit 1; }
	echo Unpacking prot.accession2taxid.gz
	gunzip -c prot.accession2taxid.gz > prot.accession2taxid
	echo Converting NR file to Kaiju database
	if [ $db_euk -eq 1 ]
	then
		gunzip -c nr.gz | convertNR -t nodes.dmp -g prot.accession2taxid -a -o kaiju_db_nr_euk.faa -l $SCRIPTDIR/taxonlist.tsv
		echo Creating BWT from Kaiju database
		mkbwt -e $exponentSA_NR -n $threadsBWT -a ACDEFGHIKLMNPQRSTVWY -o kaiju_db_nr_euk kaiju_db_nr_euk.faa
		echo Creating FM-index
		mkfmi kaiju_db_nr_euk
		echo Done!
		echo Kaiju only needs the files kaiju_db_nr_euk.fmi, nodes.dmp, and names.dmp.
		echo The remaining files can be deleted.
		echo
	else
		gunzip -c nr.gz | convertNR -t nodes.dmp -g prot.accession2taxid -a -o kaiju_db_nr.faa
		echo Creating BWT from Kaiju database
		mkbwt -e $exponentSA_NR -n $threadsBWT -a ACDEFGHIKLMNPQRSTVWY -o kaiju_db_nr kaiju_db_nr.faa
		echo Creating FM-index
		mkfmi kaiju_db_nr
		echo Done!
		echo Kaiju only needs the files kaiju_db_nr.fmi, nodes.dmp, and names.dmp.
		echo The remaining files can be deleted.
		echo
	fi
else
	echo Creating directory genomes/
	mkdir -p genomes
	if [ $db_refseq -eq 1 ]
	then
		if [ $DL -eq 1 ]
		then
			echo Downloading file list for complete genomes from RefSeq...
			wget -nv -O assembly_summary.archaea.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt
			wget -nv -O assembly_summary.bacteria.txt ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt
			awk 'BEGIN{FS="\t";OFS="/"}$12=="Complete Genome" && $11=="latest"{l=split($20,a,"/");print $20,a[l]"_genomic.gbff.gz"}' assembly_summary.bacteria.txt assembly_summary.archaea.txt > downloadlist.txt
			nfiles=`cat downloadlist.txt| wc -l`
			echo Downloading $nfiles genome files from NCBI FTP server. This may take a while...
			cat downloadlist.txt | xargs -P $parallelDL -n 1 wget -P genomes -nv
			if [ $db_viruses -eq 1 ]
			then
				echo Downloading virus genomes from RefSeq...
				wget -N -nv $wgetProgress -P genomes ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.genomic.gbff.gz
				wget -N -nv $wgetProgress -P genomes ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.genomic.gbff.gz
			fi
		fi
		if [ $db_viruses -eq 1 ]; then if [ ! -r genomes/viral.1.genomic.gbff.gz ]; then echo Missing file viral.1.genomic.gbff.gz; exit 1; fi; fi
		if [ $db_viruses -eq 1 ]; then if [ ! -r genomes/viral.2.genomic.gbff.gz ]; then echo Missing file viral.2.genomic.gbff.gz; exit 1; fi; fi

		echo Extracting protein sequences from downloaded files...
		find ./genomes -name "*.gbff.gz" | xargs -n 1 -P $parallelConversions -i gbk2faa.pl '{}' '{}'.faa
	else # must be proGenomes
		if [ $DL -eq 1 ]
		then
			echo Downloading proGenomes database...
			wget -N -nv $wgetProgress -P genomes http://progenomes.embl.de/data/repGenomes/representatives.proteins.fasta.gz
			if [ $db_viruses -eq 1 ]
			then
				echo Downloading virus genomes from RefSeq...
				wget -N -nv $wgetProgress -P genomes ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.genomic.gbff.gz
				wget -N -nv $wgetProgress -P genomes ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.genomic.gbff.gz
			fi
		fi
		if [ $db_viruses -eq 1 ]; then if [ ! -r genomes/viral.1.genomic.gbff.gz ]; then echo Missing file viral.1.genomic.gbff.gz; exit 1; fi; fi
		if [ $db_viruses -eq 1 ]; then if [ ! -r genomes/viral.2.genomic.gbff.gz ]; then echo Missing file viral.2.genomic.gbff.gz; exit 1; fi; fi

		echo Extracting protein sequences from downloaded files...
		gunzip -c genomes/representatives.proteins.fasta.gz | perl -lne 'if(/>(\d+)\.(\S+)/){print ">",$2,"_",$1}else{y/BZ/DE/;s/[^ARNDCQEGHILKMFPSTWYV]//gi;print if length}' > genomes/representatives.proteins.fasta.gz.faa
		find ./genomes -name "viral.*.gbff.gz" | xargs -n 1 -P $parallelConversions -i gbk2faa.pl '{}' '{}'.faa
	fi

	# on-the-fly substitution of taxon IDs found in merged.dmp by their updated IDs
	cat genomes/*.faa | perl -lsne 'BEGIN{open(F,$m);while(<F>){@F=split(/[\|\s]+/);$h{$F[0]}=$F[1]}}if(/(>.+)_(\d+)/){print $1,"_",defined($h{$2})?$h{$2}:$2;}else{print}' -- -m=merged.dmp  >kaiju_db.faa

	echo Creating Borrows-Wheeler transform...
	mkbwt -n $threadsBWT -e $exponentSA -a ACDEFGHIKLMNPQRSTVWY -o kaiju_db kaiju_db.faa
	echo Creating FM-Index...
	mkfmi kaiju_db
	echo Done!
	echo Kaiju only needs the files kaiju_db.fmi, nodes.dmp, and names.dmp.
	echo The remaining files and the folder genomes/ can be deleted.
fi


