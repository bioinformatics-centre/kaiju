#!/usr/bin/perl -w
# This program reads gbk files and extracts all amino acid sequences from the
# /translation fields into a FASTA file. The FASTA header contains a sequential
# number followed by the taxon id, which is extracted from the
# /db_xref="taxon:<ID>" field. Only letters in the 20 letter amino acid
# alphabet are retained in the FASTA file.
#
# author: Peter Menzel
#
# This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh
# Kaiju is licensed under the GPLv3, see the file LICENSE.
#

use strict;
use IO::Uncompress::AnyUncompress qw(anyuncompress $AnyUncompressError) ;

my $c = 0;
my $t = 0;
my $taxid;

if(!defined $ARGV[1]) { die "Usage: $0 infile.gbk outfile.faa"; }
open(OUT,">",$ARGV[1]) or die "Could not open file $ARGV[1] for writing.";
my $in_fh = new IO::Uncompress::AnyUncompress $ARGV[0] or die "Opening input file failed: $AnyUncompressError\n";

while(<$in_fh>) {
	chomp;
	if(m,/db_xref="taxon:(\d+)",) {
		$taxid = $1;
	}
	elsif(m,\s+/translation="([^"]+)",)  {  
		if(!defined($taxid)) { die "No taxon id found in gbk file $ARGV\n";}
		$c++;
		print OUT ">$c\_$taxid\n";
		my $seq = $1;
		$seq =~ tr/BZ/DE/;  # a.a. alphabet specifies `B' matches `N' or `D', and `Z' matches `Q' or `E.', here we use substitution with higher score
		$seq =~ s/[^ARNDCQEGHILKMFPSTWYV]//gi;
		print OUT "$seq\n";
	}
	elsif(m,\s+/translation="([^"]+)$,) {
		if(!defined($taxid)) { die "No taxon id found in gbk file $ARGV\n";}
		$c++;
		print OUT ">$c\_$taxid\n";
		$t = 1;
		my $seq = $1;
		$seq =~ tr/BZ/DE/;
		$seq =~ s/[^ARNDCQEGHILKMFPSTWYV]//gi;
		print OUT "$seq\n";
	}
	elsif($t) {
		if(m,",) { 
			tr/BZ/DE/;
			s/[^ARNDCQEGHILKMFPSTWYV]//gi;
			print OUT $_,"\n";   
			$t = 0;
		}
		else { 
			tr/BZ/DE/;
			s/[^ARNDCQEGHILKMFPSTWYV]//gi;
			print OUT $_,"\n";   
		}
	}
}
close($in_fh);
close(OUT);

