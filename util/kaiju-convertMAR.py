#!/usr/bin/env python
"""
This script maps downloaded sequence data from mmp_id to taxonomic lineage
using MMP metadata.

The script assumes mmp_id in column 107 and taxonomic lineage in 38.

Espen M. Robertsen, 2017. espen.m.robertsen@uit.no
"""

import argparse
import os
import sys
import json
from collections import Counter


def readfq(fp):
    """
    FAST(A/Q) generator.

    https://github.com/lh3/readfq
    """
    # this is a generator function
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break

        if not last:
            break

        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break

            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record

            if not last:
                break

        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs)

                    # yield a fastq record
                    break

            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break


def parse_tsv(tsv, db_map=None):
    """
    DEPRECATED as of 2022

    Hardcoded columns refer to taxon_lineage_ids and mmp_ID within MarRef
    and MarDB reference files. Returns a dictionary linking mmp_ID to
    last element of taxon_lineage_ids.
    """
    if db_map is None:
        db_map = dict()
    with open(tsv) as fh:
        # advance beyond header
        next(fh)
        for line in fh:
            tokens = line.strip("\r\n").split("\t")
            # column header: taxon_lineage_ids
            lineage = tokens[34].rpartition("|")[-1]
            # column header: mmp_ID
            db_map[tokens[101]] = lineage
    return db_map

def parse_json(metadata, db_map=None):
    """
    Parse json files retrieved from mmp api. Returns a dictionary linking NCBI
    assembly accession (GCA) to NCBI taxonomic id
    """
    if db_map is None:
        tmp = json.load(open(metadata, "r"))
        db_map = {entry['x']:entry['y_yName'][0].split(":")[-1] for entry in tmp['graph']}
    else:
        tmp = json.load(open(metadata, "r"))
        tmp = {entry['x']:entry['y_yName'][0].split(":")[-1] for entry in tmp['graph']}
        db_map = {**db_map, **tmp}
    return db_map

def parse_nodes(nodes):
    """
    Parse NCBI taxonomy's nodes.dmp, grabbing the taxonomy identifier in the
    first column. Returns a set containing the IDs.
    """
    ids = set()
    with open(nodes) as fh:
    	for line in fh:
    		ids.add(line.partition("|")[0].strip())
    return ids

def process_genomes(genome_dir, lineage_map, tax_ids):
    """
    Walk through the genomes directory parsing .faa files, updating the
    sequence header to include the NCBI taxonomy. Prints to STDOUT.
    """
    processed_accessions = set()
    warnings = Counter()
    for root, dirs, files in os.walk(genome_dir):
        for dirname in dirs:
            accession = dirname
            filename = os.path.join(dirname, "protein.faa")
            filename = os.path.join(root, filename)
            with open(filename) as fh:
                for i, (name, seq, _) in enumerate(readfq(fh), start=1):

                    # If this mmp accession is a duplicate (already processed), skip it.
                    if accession in processed_accessions:
                        warnings.update(["duplicate"])
                        break

                    # If this mmp accession has no taxonomic linage for some reason, skip it.
                    if not accession in lineage_map:
                        warnings.update(["nolineage"])
                        break

                    # If taxonomic id is not in nodes.dmp, skip.
                    if not lineage_map[accession] in tax_ids:
                        warnings.update(["notax"])
                        break

                    # If sequence contains dubious characters, skip.
                    if '*' in seq:
                        warnings.update(["asterix"])
                        continue

                    current_lineage = lineage_map[accession]
                    sys.stdout.write(
                        ">{index}_{id}_{lineage}\n{seq}\n".format(
                            index=i, id=accession, lineage=current_lineage, seq=seq
                        )
                    )
                processed_accessions.add(accession)
    return warnings

def process_mardb(ref, nodes, genomes):
    """
    Parse TSVs and reformat MarDB genomes to accommodate Kaiju FASTA
    reference requirements.
    """
    if args.ref:
        lineage_map = parse_json(ref)
    else:
        exit("No JSON metadata specified (--ref)")
    #lineage_map = parse_json(db, db_map=lineage_map)
    tax_ids = parse_nodes(nodes)
    warn = process_genomes(genomes, lineage_map, tax_ids)
    # Output simple statistic to stderr
    sys.stderr.write(
        (
            "# Warnings / Sequences removed:\n"
            "# Duplicates: {duplicates}\n"
            "# Accessions with no taxonomic lineage: {no_lineage}\n"
            "# Accessions with no taxid in nodes.dmp: {no_tax}\n"
            "# Sequences with asterix: {asterix}\n"
        ).format(
            duplicates=warn["duplicate"],
            no_lineage=warn["nolineage"],
            no_tax=warn["notax"],
            asterix=warn["asterix"],
        )
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--ref", default=False, help="Mar metadta JSON file path")
    #p.add_argument("--db", default=False, help="MarDB TSV file path")
    p.add_argument("--nodes", default="nodes.dmp", help="NCBI nodes.dmp file path")
    p.add_argument("--genomes", default="genomes", help="genomes download directory")
    args = p.parse_args()
    process_mardb(args.ref, args.nodes, args.genomes)
