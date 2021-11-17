#!/bin/bash

set -e

echo "Fetching test data"
git clone https://github.com/pmenzel/kaiju-testdata.git

echo "Running kaiju"
./bin/kaiju -z 2 -t nodes.dmp -f viruses/kaiju_db_viruses.fmi -i kaiju-testdata/sars-cov-2_1.fastq.gz -o sars-cov-2_1.out
./bin/kaiju -z 2 -t nodes.dmp -f viruses/kaiju_db_viruses.fmi -i kaiju-testdata/sars-cov-2_2.fastq.gz -o sars-cov-2_2.out

echo "Running kaiju-multi"
./bin/kaiju-multi -z 2 -t nodes.dmp -f viruses/kaiju_db_viruses.fmi -i kaiju-testdata/sars-cov-2_1.fastq.gz,kaiju-testdata/sars-cov-2_2.fastq.gz > multi-sars-cov-2_combined.out
./bin/kaiju-multi -z 2 -t nodes.dmp -f viruses/kaiju_db_viruses.fmi -i kaiju-testdata/sars-cov-2_1.fastq.gz,kaiju-testdata/sars-cov-2_2.fastq.gz -o multi-sars-cov-2_1.out,multi-sars-cov-2_2.out

echo "Running kaiju2table"
./bin/kaiju2table -t nodes.dmp -n names.dmp -e -r species -o sars-cov-2_1.table sars-cov-2_1.out
./bin/kaiju2table -t nodes.dmp -n names.dmp -e -r species -o sars-cov-2_2.table sars-cov-2_2.out
./bin/kaiju2table -t nodes.dmp -n names.dmp -e -r species -o sars-cov-2_combined.table sars-cov-2_1.out sars-cov-2_2.out
./bin/kaiju2table -t nodes.dmp -n names.dmp -e -r species -o multi-sars-cov-2_1.table multi-sars-cov-2_1.out
./bin/kaiju2table -t nodes.dmp -n names.dmp -e -r species -o multi-sars-cov-2_2.table multi-sars-cov-2_2.out
./bin/kaiju2table -t nodes.dmp -n names.dmp -e -r species -o multi-sars-cov-2_combined.table multi-sars-cov-2_combined.out

echo "Testing output files"
# kaiju-multi should be same output in individual output files and a combined outputfile from stdout
cmp <(cat multi-sars-cov-2_combined.out | sort) <(cat multi-sars-cov-2_1.out multi-sars-cov-2_2.out | sort)
# the sum of read counts per taxon from single tables should equal the counts in the table from combined output
cmp <(perl -lsane '$h{$F[3]}+=$F[2]}{map { print "$_ $h{$_}";} sort keys %h' multi-sars-cov-2_1.table multi-sars-cov-2_2.table) <(perl -lsane '$h{$F[3]}+=$F[2]}{map { print "$_ $h{$_}";} sort keys %h' multi-sars-cov-2_combined.table)
# tables from single samples using normal kaiju should be the same as from kaiju-multi:
cmp <(sed '/^multi-/s/^multi-//' multi-sars-cov-2_1.table) sars-cov-2_1.table
cmp <(sed '/^multi-/s/^multi-//' multi-sars-cov-2_2.table) sars-cov-2_2.table

