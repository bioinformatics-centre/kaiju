#!/usr/bin/python

# This script maps downloaded sequence data from mmp_id to taxonomic lineage using MMP metadata. 
# The script assumes mmp_id in column 107 and taxonomic lineage in 38. 

from collections import defaultdict
import os, sys, HTSeq

mapdict = defaultdict(str)

# Map downloaded metadata to the defaultdict 'mapdict'
with open ('MarRef.tsv') as marref:
	for line in marref:
		if line.startswith("base_ID"):
			continue
		records = line.strip().split("\t")
		lineage = records[37].split('|')
		lineage = lineage[-1]
		mapdict[records[106]] = lineage

with open ('MarDB.tsv') as mardb:
	for line in mardb:
		if line.startswith("base_ID"):
			continue
		records = line.strip().split("\t")
		lineage = records[37].split('|')
		lineage = lineage[-1]
		mapdict[records[106]] = lineage

# Verify taxids against nodes.dmp
taxids = set()
with open ('nodes.dmp') as nodes:
	for line in nodes:
		data = line.split("|")
		taxids.add(data[0].strip())

# Map mmp_id from downloaded sequences to lineages and write to stdout.
processed_mmp = set()
warnings = defaultdict(int)
for root, dirs, files in os.walk("./genomes"):
	for filename in files:
		if filename.endswith(".faa"):
			counter = 1
			for seq in HTSeq.FastaReader(os.path.join(root,filename)):
				mmp_id = seq.name.split("_")
				mmp_id = mmp_id[-1]

				# If this mmp accession is a duplicate (already processed), skip it.
				if mmp_id in processed_mmp:
					warnings['duplicate']+=1
					break
				# If this mmp accession has no taxonomic linage for some reason, skip it.
				if not mapdict[mmp_id]:
					warnings['nolineage']+=1
					break
				# If taxonomic id is not in nodes.dmp, skip.
				if not mapdict[mmp_id] in taxids:
					warnings['notax']+=1
					break
				# If sequence contains dubious characters, skip.
				if '*' in seq.seq:
					warnings['asterix']+=1
					continue 

				lineage = mapdict[mmp_id]
				header = str(counter)+"_"+mmp_id+"_"+str(lineage)
				counter += 1
				print ">"+header+"\n", seq.seq
			processed_mmp.add(mmp_id)

sys.stderr.write("Warnings / Sequences removed:\nDuplicates: "+str(warnings['duplicate'])+"\nAccessions with no taxonomic lineage: "+str(warnings['nolineage'])+"\nAccessions with no taxid in nodes.dmp: "+str(warnings['notax'])+"\nSequences with asterix: "+str(warnings['asterix'])+"\n") 
