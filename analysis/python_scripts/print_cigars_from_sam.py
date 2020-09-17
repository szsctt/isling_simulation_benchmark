#!/usr/bin/env python3

# this script prints cigar strings from reads
# usage: python3 print_reads_from_sam.py <samfile> <read_id_1> ...

# if read id ends in /1 or /2, print only read1 or read2 matching that qname, otherwise print both

from sys import argv
import pysam
import pdb

def main(argv):
	# check argv length
	if len(argv)  < 3:
		print("usage: python3 print_reads_from_sam.py <samfile> <read_id_1> <read_id_2> ...")
		return()
	
	# build dict of reads we're looking for
	reads = {"1": [], "2":[]}
	for read in argv[2:]:
		if read[-2:] == "/1":
			reads["1"].append(read[:-2])
		elif read[-2:] == "/2":
			reads["2"].append(read[:-2])
		else:
			reads["1"].append(read)
			reads["2"].append(read)
	
	# loop through samfile and print cigars for matches
	found = []
	samfile = pysam.AlignmentFile(argv[1])
	for read in samfile:
		if read.is_reverse:
			ori = "reverse"
		else:
			ori = "forward"
		if read.is_read1 is True:
			if read.qname in reads["1"]:
				if read.is_secondary or read.is_supplementary:
					continue
				else:
					print(f"{read.qname}/1 ({ori}): {read.cigarstring}")
					found.append(f"{read.qname}/1")
		elif read.is_read2 is True:
			if read.qname in reads["2"]:
				if read.is_secondary or read.is_supplementary:
					continue
				else:
					print(f"{read.qname}/2 ({ori}): {read.cigarstring}")
					found.append(f"{read.qname}/2")
	
	# check if there are any reads we didn't find
	for read_num in reads:
		for readid in reads[read_num]:
			if f"{readid}/{read_num}" not in found:
				print(f"didn't find read {readid}/{read_num} in sam file {argv[1]}")

	
	
if __name__ == "__main__":
	main(argv)	


