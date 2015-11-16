#!/usr/bin/perl -w
# This program reads gbk files and extracts all amino acid sequences from the
# /translation fields into a FASTA file. The FASTA header contains a sequential
# number followed by the taxon id, which is extracted from the
# /db_xref="taxon:<ID>" field. Only letters in the 20 letter amino acid
# alphabet are retained in the FASTA file.
#
# author: Peter Menzel
#
# This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh
# Kaiju is licensed under the GPLv3, see the file LICENSE.
#

use strict;
my $c = 0;
my $t = 0;
my $taxid;

if(!defined $ARGV[1]) { die "Usage: $0 infile.gbk outfile.faa"; }
open(IN,$ARGV[0]) or die "Could not open file $ARGV[0] for reading.";
open(OUT,">",$ARGV[1]) or die "Could not open file $ARGV[1] for writing.";

while(<IN>) {
	chomp;
	if(m,/db_xref="taxon:(\d+)",) {
		$taxid = $1;
	}
	elsif(m,\s+/translation="([^"]+)",)  {  
		if(!defined($taxid)) { die "No taxon id found in gbk file $ARGV\n";}
		$c++;
		print OUT ">$c\_$taxid\n";
		my $seq = $1;
		$seq =~ s/[^ARNDCQEGHILKMFPSTWYV]//gi;
		print OUT "$seq\n";
	}
	elsif(m,\s+/translation="([^"]+)$,) {
		if(!defined($taxid)) { die "No taxon id found in gbk file $ARGV\n";}
		$c++;
		print OUT ">$c\_$taxid\n";
		$t = 1;
		my $seq = $1;
		$seq =~ s/[^ARNDCQEGHILKMFPSTWYV]//gi;
		print OUT "$seq\n";
	}
	elsif($t) {
		if(m,",) { 
			s/[^ARNDCQEGHILKMFPSTWYV]//gi;
			print OUT $_,"\n";   
			$t = 0;
		}
		else { 
			s/[^ARNDCQEGHILKMFPSTWYV]//gi;
			print OUT $_,"\n";   
		}
	}
}
close(IN);
close(OUT);
exit;

