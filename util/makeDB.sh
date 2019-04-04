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
db_mar_ref=0
db_mar_db=0
db_mar_mags=0
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
echo  "$s" -m  Marine Metagenomics Portal \(MMP\) marine reference databases, MarRef, MarDB
echo  "$tab"   and MarDB Metagenomic assembled genomes \(MAGS\) \(https://mmp.sfb.uit.no\)
echo  "$tab"   To build discrete divisions, provide -m \<division\>, where divisions can be
echo  "$tab"   \'MarRef\', \'MarDB\' or \'MarDBMAGS\'. No arguments will build all divisions.
echo  "$tab"   Multiple divisions separated by space can be selected.
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
			# If no arguments to -m, or next element in $@ is a d flag, build all 3 divisions
			if [ -z "$2" ] || case $2 in -*) ;; *) false;; esac; then
				db_mar=1
				db_mar_ref=1
				db_mar_db=1
				db_mar_mags=1
				# If arguments are supplied, loop through the rest of $@ and set individual discrete divisions
			else
				db_mar=1
				for param in $@; do
					if [ $param = 'MarRef' ]; then
						db_mar_ref=1
						shift
						continue
					fi
					if [ $param = 'MarDB' ]; then
						db_mar_db=1
						shift
						continue
					fi
					if [ $param = 'MarDBMAGS' ]; then
						db_mar_mags=1
						shift
						continue
					fi
				done
			fi
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

[ $db_plasmids -eq 1 -o $db_viruses -eq 1 -o $db_refseq -eq 1 -o $db_progenomes -eq 1 -o $db_nr -eq 1 -o $db_euk -eq 1 -o $db_mar -eq 1 ] || { echo "Error: Use one of the options -r, -p, -n, -v, -l, -m, or -e"; usage; exit 1; }

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
command -v kaiju-gbk2faa.pl >/dev/null 2>/dev/null || { echo Error: kaiju-gbk2faa.pl not found in $PATH; exit 1; }
command -v kaiju-mkfmi >/dev/null 2>/dev/null || { echo Error: kaiju-mkfmi not found in $PATH; exit 1; }
command -v kaiju-mkbwt >/dev/null 2>/dev/null || { echo Error: kaiju-mkbwt not found in $PATH; exit 1; }
command -v kaiju-convertNR >/dev/null 2>/dev/null || { echo Error: kaiju-convertNR not found in $PATH; exit 1; }

if [ $db_euk -eq 1 ]
then
	[ -r $SCRIPTDIR/kaiju-taxonlistEuk.tsv ] || { echo Error: file kaiju-taxonlistEuk.tsv not found in $SCRIPTDIR; exit 1; }
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
		echo Creating directory genomes/
                mkdir -p genomes
		echo Downloading MarRef reference genomes from the Marine Metagenomics Portal \(MMP\). This may take a while...
		if [ $db_mar_ref -eq 1 ]; then
		    wget -nv -O dl_list_marref_protein.txt https://s1.sfb.uit.no/public/mar/Resources/kaiju/dl_list_marref_protein.txt
		    cat dl_list_marref_protein.txt | xargs -P $parallelDL wget -P genomes -q || true
		fi
		echo Downloading MarDB complete genomes from the Marine Metagenomics Portal \(MMP\). This may take a while...
		if [ $db_mar_db -eq 1 ]; then
		    wget -nv -O dl_list_mardb_no_mags_protein.txt https://s1.sfb.uit.no/public/mar/Resources/kaiju/dl_list_mardb_no_mags_protein.txt
		    cat dl_list_mardb_no_mags_protein.txt | xargs -P $parallelDL wget -P genomes -q || true
		fi
		echo Downloading MarDBMAGS metagenomic assembled genomes from the Marine Metagenomics Portal \(MMP\). This may take a while...
		if [ $db_mar_mags -eq 1 ]; then
		    wget -nv -O dl_list_mardb_mags_protein.txt https://s1.sfb.uit.no/public/mar/Resources/kaiju/dl_list_mardb_mags_protein.txt
		    cat dl_list_mardb_mags_protein.txt | xargs -P $parallelDL wget -P genomes -q || true
		fi
		echo Downloading necessary metadata from MMP
		wget -nv -O MarRef.tsv https://s1.sfb.uit.no/public/mar/MarRef/Metadatabase/Current.tsv
		wget -nv -O MarDB.tsv https://s1.sfb.uit.no/public/mar/MarDB/Metadatabase/Current.tsv
	fi
	[ -r MarRef.tsv ] || { echo Missing file MarRef.tsv; exit 1; }
	[ -r MarDB.tsv ] || { echo Missing file MarDB.tsv; exit 1; }
	[ -r $SCRIPTDIR/kaiju-convertMAR.py ] || { echo Error: file kaiju-convertMAR.py not found in $SCRIPTDIR; exit 1; }
	echo Converting MMP data to kaiju format
	python $SCRIPTDIR/kaiju-convertMAR.py > kaiju_db_tmp.faa
	echo On the fly substitution with merged.dmp
	cat kaiju_db_tmp.faa | perl -lsne 'BEGIN{open(F,$m);while(<F>){@F=split(/[\|\s]+/);$h{$F[0]}=$F[1]}}if(/(>.+)_(\d+)/){print $1,"_",defined($h{$2})?$h{$2}:$2;}else{print}' -- -m=merged.dmp > kaiju_db.faa
	rm kaiju_db_tmp.faa
	echo Building Kaiju reference
	kaiju-mkbwt -n $threadsBWT -e $exponentSA -a ACDEFGHIKLMNPQRSTVWY -o kaiju_db kaiju_db.faa
	kaiju-mkfmi kaiju_db
	exit 0
fi


if [ $db_nr -eq 1 -o $db_euk -eq 1 ]
then
	if [ $DL -eq 1 ]
	then
		echo Downloading file nr.gz
		wget -c -N -nv $wgetProgress ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
		echo Downloading file prot.accession2taxid.gz
		wget -c -N -nv $wgetProgress ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
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
			gunzip -c nr.gz | kaiju-convertNR -t nodes.dmp -g prot.accession2taxid -a -o kaiju_db_nr_euk.faa -l $SCRIPTDIR/kaiju-taxonlistEuk.tsv
		fi
		[ -r kaiju_db_nr_euk.faa ] || { echo Missing file kaiju_db_nr_euk.faa; exit 1; }
		echo Creating BWT from Kaiju database
		kaiju-mkbwt -e $exponentSA_NR -n $threadsBWT -a ACDEFGHIKLMNPQRSTVWY -o kaiju_db_nr_euk kaiju_db_nr_euk.faa
		echo Creating FM-index
		kaiju-mkfmi kaiju_db_nr_euk
		echo Done!
		echo Kaiju only needs the files kaiju_db_nr_euk.fmi, nodes.dmp, and names.dmp.
		echo The remaining files can be deleted.
		echo
	else
		if [ $index_only -eq 0 ]
		then
			echo Converting NR file to Kaiju database
			gunzip -c nr.gz | kaiju-convertNR -t nodes.dmp -g prot.accession2taxid -a -o kaiju_db_nr.faa
		fi
		[ -r kaiju_db_nr.faa ] || { echo Missing file kaiju_db_nr.faa; exit 1; }
		echo Creating BWT from Kaiju database
		kaiju-mkbwt -e $exponentSA_NR -n $threadsBWT -a ACDEFGHIKLMNPQRSTVWY -o kaiju_db_nr kaiju_db_nr.faa
		echo Creating FM-index
		kaiju-mkfmi kaiju_db_nr
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
			find ./genomes -name "*.gbff.gz" | xargs -n 1 -P $parallelConversions -IXX kaiju-gbk2faa.pl XX XX.faa
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
			find ./genomes -name "viral.*.gbff.gz" | xargs -n 1 -P $parallelConversions -IXX kaiju-gbk2faa.pl XX XX.faa
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
			find ./genomes -name "viral.*.gbff.gz" | xargs -n 1 -P $parallelConversions -IXX kaiju-gbk2faa.pl XX XX.faa
		elif [ $db_plasmids -eq 1 ]
		then
			if [ $DL -eq 1 ]
			then
				echo Downloading plasmid genomes from RefSeq...
				wget -N -nv $wgetProgress -P genomes ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.1.genomic.gbff.gz
				wget -N -nv $wgetProgress -P genomes ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.2.genomic.gbff.gz
				wget -N -nv $wgetProgress -P genomes ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.3.genomic.gbff.gz
				wget -N -nv $wgetProgress -P genomes ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.4.genomic.gbff.gz
				wget -N -nv $wgetProgress -P genomes ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plasmid/plasmid.5.genomic.gbff.gz
			fi
			if [ ! -r genomes/plasmid.1.genomic.gbff.gz ]; then echo Missing file plasmid.1.genomic.gbff.gz; exit 1; fi;
			if [ ! -r genomes/plasmid.2.genomic.gbff.gz ]; then echo Missing file plasmid.2.genomic.gbff.gz; exit 1; fi;
			if [ ! -r genomes/plasmid.3.genomic.gbff.gz ]; then echo Missing file plasmid.3.genomic.gbff.gz; exit 1; fi;
			if [ ! -r genomes/plasmid.4.genomic.gbff.gz ]; then echo Missing file plasmid.4.genomic.gbff.gz; exit 1; fi;
			if [ ! -r genomes/plasmid.5.genomic.gbff.gz ]; then echo Missing file plasmid.5.genomic.gbff.gz; exit 1; fi;
			echo Extracting protein sequences from downloaded files...
			find ./genomes -name "plasmid.*.gbff.gz" | xargs -n 1 -P $parallelConversions -IXX kaiju-gbk2faa.pl XX XX.faa
		fi

		# on-the-fly substitution of taxon IDs found in merged.dmp by their updated IDs
		#old, too many arguments, see #89: cat genomes/*.faa | perl -lsne 'BEGIN{open(F,$m);while(<F>){@F=split(/[\|\s]+/);$h{$F[0]}=$F[1]}}if(/(>.+)_(\d+)/){print $1,"_",defined($h{$2})?$h{$2}:$2;}else{print}' -- -m=merged.dmp  >kaiju_db.faa
		#new, but untested on Mac:
		find genomes -name '*.faa' -print0 | xargs -0 cat |  perl -lsne 'BEGIN{open(F,$m);while(<F>){@F=split(/[\|\s]+/);$h{$F[0]}=$F[1]}}if(/(>.+)_(\d+)/){print $1,"_",defined($h{$2})?$h{$2}:$2;}else{print}' -- -m=merged.dmp  >kaiju_db.faa
	fi

	[ -r kaiju_db.faa ] || { echo Missing file kaiju_db.faa; exit 1; }
	echo Creating Borrows-Wheeler transform...
	kaiju-mkbwt -n $threadsBWT -e $exponentSA -a ACDEFGHIKLMNPQRSTVWY -o kaiju_db kaiju_db.faa
	echo Creating FM-Index...
	kaiju-mkfmi kaiju_db
	echo Done!
	echo Kaiju only needs the files kaiju_db.fmi, nodes.dmp, and names.dmp.
	echo The remaining files and the folder genomes/ can be deleted.
fi
