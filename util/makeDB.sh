#!/bin/sh
#
# This file is part of Kaiju, Copyright 2015-2018 Peter Menzel and Anders Krogh
# Kaiju is licensed under the GPLv3, see the file LICENSE.
#
SCRIPTDIR=$(dirname $0)

PATH=$SCRIPTDIR:$PATH

db_viruses=0
db_refseq=0
db_progenomes=0
db_nr=0
db_euk=0
db_mar=0
db_plasmids=0
threadsBWT=5
parallelDL=5
parallelConversions=5
exponentSA=3
exponentSA_NR=5
DL=1
index_only=0

usage() {
s=" "
tab="    "
echo
echo This program creates a protein reference database and index for Kaiju.
echo Several source databases can be used and one of these options must be set:
echo
echo  "$s" -r  All complete bacterial and archaeal genomes in the NCBI RefSeq database
echo
echo  "$s" -p  All proteins belonging to the set of representative genomes
echo  "$tab"   from the proGenomes database
echo
echo  "$s" -n  NCBI BLAST non-redundant protein database \"nr\":
echo  "$tab"   only Archaea, bacteria, and viruses
echo
echo  "$s" -e  NCBI BLAST non-redundant protein database \"nr\":
echo  "$tab"   like -n, but additionally including fungi and microbial eukaryotes
echo
echo  "$s" -m  Marine Metagenomics Portal \(MMP\) marine reference databases, MarRef and MarDB
echo  "$tab"   \(https://mmp.sfb.uit.no\)
echo
echo  "$s" -v   Viral genomes from RefSeq, can also be used together with -n or -p
echo
echo  "$s" -l   Plasmid genomes from RefSeq
echo
echo Additional options:
echo
echo  "$s" -t X  Set number of parallel threads for index construction to X \(default:5\)
echo  "       The more threads are used, the higher the memory requirement becomes."
echo
echo  "$s" --noDL  Do not download files, but use the existing files in the folder.
echo
echo  "$s" --index-only  Only create BWT and FMI from kaiju_db*.faa files, implies --noDL.
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
        --index-only)
            index_only=1
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
        -l|--plasmids)
            db_plasmids=1
            ;;
        -p|--progenomes)
            db_progenomes=1
            ;;
        -r|--refseq)
            db_refseq=1
            ;;
        -m|--mardb)
            db_mar=1
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

[ $db_plasmids -eq 1 -o $db_viruses -eq 1 -o $db_refseq -eq 1 -o $db_progenomes -eq 1 -o $db_nr -eq 1 -o $db_euk -eq 1 -o $db_mar -eq 1 ] || { echo "Error: Use one of the options -r, -p, -n, -v, -m or -e"; usage; exit 1; }

#check if necessary programs are in the PATH
command -v awk >/dev/null 2>/dev/null || { echo Error: awk not found; exit 1; }
command -v wget >/dev/null 2>/dev/null || { echo Error: wget not found; exit 1; }
command -v xargs >/dev/null 2>/dev/null || { echo Error: xargs not found; exit 1; }
command -v tar >/dev/null 2>/dev/null || { echo Error: tar not found; exit 1; }
command -v gunzip >/dev/null 2>/dev/null || { echo Error: gunzip not found; exit 1; }
command -v perl >/dev/null 2>/dev/null || { echo Error: perl not found; exit 1; }
command -v python >/dev/null 2>/dev/null || { echo Error: python not found; exit 1; }

#test if option --show-progress is available for wget, then use it when downloading
wgetProgress=""
wget --help | grep -q -- --show-progress && wgetProgress='--show-progress'

#check that all programs from Kaiju are usable
command -v gbk2faa.pl >/dev/null 2>/dev/null || { echo Error: gbk2faa.pl not found in $PATH; exit 1; }
command -v mkfmi >/dev/null 2>/dev/null || { echo Error: mkfmi not found in $PATH; exit 1; }
command -v mkbwt >/dev/null 2>/dev/null || { echo Error: mkbwt not found in $PATH; exit 1; }
command -v convertNR >/dev/null 2>/dev/null || { echo Error: convertNR not found in $PATH; exit 1; }
command -v convert_mar_to_kaiju.py >/dev/null 2>/dev/null || { echo Error: convert_mar_to_kaiju.py not found in $PATH; exit 1; }

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
[ -r taxdump.tar.gz ] || { echo Missing file taxdump.tar.gz; exit 1; }
echo Extracting file taxdump.tar.gz
tar xf taxdump.tar.gz nodes.dmp names.dmp merged.dmp

if [ $db_mar -eq 1 ]
then
	if [ $DL -eq 1 ]
	then
		echo Downloading list of marine genomes from the Marine Metagenomics Portal \(MMP\)
		wget -nv -O download_list.txt https://s1.sfb.uit.no/public/mar/Resources/kaiju/download_list.txt
		echo Downloading necessary metadata from MMP
		wget -nv -O MarRef.tsv https://s1.sfb.uit.no/public/mar/MarRef/Metadatabase/Current.tsv
		wget -nv -O MarDB.tsv https://s1.sfb.uit.no/public/mar/MarDB/Metadatabase/Current.tsv
		echo Creating directory genomes/
		mkdir -p genomes
		echo Downloading Mar reference genomes from the MMP. This may take a while...
		cat download_list.txt | xargs -P $parallelDL wget -P genomes -nv
	fi
	[ -r download_list.txt ] || { echo Missing file download_list.txt; exit 1; }
	[ -r MarRef.tsv ] || { echo Missing file MarRef.tsv; exit 1; }
	[ -r MarDB.tsv ] || { echo Missing file MarDB.tsv; exit 1; }
	[ -r $SCRIPTDIR/convert_mar_to_kaiju.py ] || { echo Error: file convert_mar_to_kaiju.py not found in $SCRIPTDIR; exit 1; }
	echo Converting MMP data to kaiju format
	python $SCRIPTDIR/convert_mar_to_kaiju.py > kaiju_db_tmp.faa
	echo On the fly substitution with merged.dmp
	cat kaiju_db_tmp.faa | perl -lsne 'BEGIN{open(F,$m);while(<F>){@F=split(/[\|\s]+/);$h{$F[0]}=$F[1]}}if(/(>.+)_(\d+)/){print $1,"_",defined($h{$2})?$h{$2}:$2;}else{print}' -- -m=merged.dmp > kaiju_db.faa
	rm kaiju_db_tmp.faa
	echo Building Kaiju reference
	mkbwt -n $threadsBWT -e $exponentSA -a ACDEFGHIKLMNPQRSTVWY -o kaiju_db kaiju_db.faa
	mkfmi kaiju_db
	exit 0
fi


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
	if [ $db_euk -eq 1 ]
	then
		if [ $index_only -eq 0 ]
		then
			echo Converting NR file to Kaiju database
			gunzip -c nr.gz | convertNR -t nodes.dmp -g prot.accession2taxid -a -o kaiju_db_nr_euk.faa -l $SCRIPTDIR/taxonlist.tsv
		fi
		[ -r kaiju_db_nr_euk.faa ] || { echo Missing file kaiju_db_nr_euk.faa; exit 1; }
		echo Creating BWT from Kaiju database
		mkbwt -e $exponentSA_NR -n $threadsBWT -a ACDEFGHIKLMNPQRSTVWY -o kaiju_db_nr_euk kaiju_db_nr_euk.faa
		echo Creating FM-index
		mkfmi kaiju_db_nr_euk
		echo Done!
		echo Kaiju only needs the files kaiju_db_nr_euk.fmi, nodes.dmp, and names.dmp.
		echo The remaining files can be deleted.
		echo
	else
		if [ $index_only -eq 0 ]
		then
			echo Converting NR file to Kaiju database
			gunzip -c nr.gz | convertNR -t nodes.dmp -g prot.accession2taxid -a -o kaiju_db_nr.faa
		fi
		[ -r kaiju_db_nr.faa ] || { echo Missing file kaiju_db_nr.faa; exit 1; }
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
	if [ $index_only -eq 0 ]
	then
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
			find ./genomes -name "*.gbff.gz" | xargs -n 1 -P $parallelConversions -IXX gbk2faa.pl XX XX.faa
		elif [ $db_progenomes -eq 1 ]
		then
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
			find ./genomes -name "viral.*.gbff.gz" | xargs -n 1 -P $parallelConversions -IXX gbk2faa.pl XX XX.faa
		elif [ $db_viruses -eq 1 ]
		then
			if [ $DL -eq 1 ]
			then
				echo Downloading virus genomes from RefSeq...
				wget -N -nv $wgetProgress -P genomes ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.genomic.gbff.gz
				wget -N -nv $wgetProgress -P genomes ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.2.genomic.gbff.gz
			fi
			if [ ! -r genomes/viral.1.genomic.gbff.gz ]; then echo Missing file viral.1.genomic.gbff.gz; exit 1; fi;
			if [ ! -r genomes/viral.2.genomic.gbff.gz ]; then echo Missing file viral.2.genomic.gbff.gz; exit 1; fi;
			echo Extracting protein sequences from downloaded files...
			find ./genomes -name "viral.*.gbff.gz" | xargs -n 1 -P $parallelConversions -IXX gbk2faa.pl XX XX.faa
		elif [ $db_plasmids -eq 1 ]
		then
			if [ $DL -eq 1 ]
			then
				echo Downloading plasmid genomes from RefSeq...
				wget -N -nv $wgetProgress -P genomes ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.1.genomic.gbff.gz
				wget -N -nv $wgetProgress -P genomes ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.2.genomic.gbff.gz
				wget -N -nv $wgetProgress -P genomes ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.3.genomic.gbff.gz
				wget -N -nv $wgetProgress -P genomes ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.4.genomic.gbff.gz
			fi
			if [ ! -r genomes/plasmid.1.genomic.gbff.gz ]; then echo Missing file plasmid.1.genomic.gbff.gz; exit 1; fi;
			if [ ! -r genomes/plasmid.2.genomic.gbff.gz ]; then echo Missing file plasmid.2.genomic.gbff.gz; exit 1; fi;
			if [ ! -r genomes/plasmid.3.genomic.gbff.gz ]; then echo Missing file plasmid.3.genomic.gbff.gz; exit 1; fi;
			if [ ! -r genomes/plasmid.4.genomic.gbff.gz ]; then echo Missing file plasmid.4.genomic.gbff.gz; exit 1; fi;
			echo Extracting protein sequences from downloaded files...
			find ./genomes -name "plasmid.*.gbff.gz" | xargs -n 1 -P $parallelConversions -IXX gbk2faa.pl XX XX.faa
		fi

		# on-the-fly substitution of taxon IDs found in merged.dmp by their updated IDs
		cat genomes/*.faa | perl -lsne 'BEGIN{open(F,$m);while(<F>){@F=split(/[\|\s]+/);$h{$F[0]}=$F[1]}}if(/(>.+)_(\d+)/){print $1,"_",defined($h{$2})?$h{$2}:$2;}else{print}' -- -m=merged.dmp  >kaiju_db.faa
	fi

	[ -r kaiju_db.faa ] || { echo Missing file kaiju_db.faa; exit 1; }
	echo Creating Borrows-Wheeler transform...
	mkbwt -n $threadsBWT -e $exponentSA -a ACDEFGHIKLMNPQRSTVWY -o kaiju_db kaiju_db.faa
	echo Creating FM-Index...
	mkfmi kaiju_db
	echo Done!
	echo Kaiju only needs the files kaiju_db.fmi, nodes.dmp, and names.dmp.
	echo The remaining files and the folder genomes/ can be deleted.
fi
