#!/usr/bin/env python3

# usage: python3 make_isling_config.py <sim_config> <output_config>


# given simulation config file, make config file for running isling on one sample from one experiment
# config file for isling looks like this:

#exp:
#  read_folder: "test/reads/"
#  out_dir: "out/pipeline-test"
#  samples:
#  	- "sample1"
#  R1_suffix: "1.fq"
#  R2_suffix: "2.fq"
#  read1-adapt: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
#  read2-adapt: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
#  mean-frag-len: "estimate"
#  merge: False
#  trim: False
#  host_name: "host"
#  host_fasta: "test/references/test_human.fa"
#  virus_name: "rep68"
#  virus_fasta: "test/references/test_AAV.fa"
#  dedup: False
#  clip-cutoff: 20
#  min-mapq: 10
#  cigar-tol: 3
#  post:
#    - filter
#    - dedup
#  merge-dist: 1000
#  min-n-merge: 1

# assume that all experiments have the same host and virus

import sys
import yaml
import pdb
import os

def main(args):

	# get inputs
	config, out = args[1:]
	
	# import simulation config yaml
	with open(config, 'r') as stream:
		try:
			in_config = yaml.safe_load(stream)
		except yaml.YAMLError as exc:
			print(exc)	
	
	# apply global options to in config
	if 'global' in in_config:		
		# get default (global) options
		default = in_config.pop('global')
		for dataset in in_config:
			for key in default:
				if key not in in_config[dataset]:
					in_config[dataset][key] = default[key]
	
	exp = list(in_config.keys())[0]
	
	# construct output yaml for isling
	out_config = {'snakedir': '.', exp: {}}
	
	# things that are always the same for these experiments
	out_config[exp]['merge'] = False
	out_config[exp]['trim'] = False
	out_config[exp]['dedup'] = False	
	out_config[exp]['clip-cutoff'] = 20
	out_config[exp]['min-mapq'] = 10
	out_config[exp]['cigar-tol'] = 3
	out_config[exp]['post'] = ['filter']
	out_config[exp]['merge-dist'] = 1000
	out_config[exp]['merge-n-min'] = 1
	out_config[exp]['R1_suffix'] = '1.fq'
	out_config[exp]['R2_suffix'] = '2.fq'
	out_config[exp]['read1-adapt'] = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
	out_config[exp]['read2-adapt'] = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
	
	# get read folder
	out_config[exp]['read_folder'] = os.path.join(in_config[exp]['out_directory'], exp, 'sim_reads')
	
	# set output directory
	out_config[exp]['out_dir'] = os.path.join(in_config[exp]['out_directory'])	
	
	# set sample
	out_config[exp]['samples'] = ['cond0.rep0']
	
	# host
	host = list(in_config[exp]['hosts'].keys())[0]
	out_config[exp]['host_name'] = host
	out_config[exp]['host_fasta'] = in_config[exp]['hosts'][host]
	
	# virus
	virus = list(in_config[exp]['viruses'].keys())[0]
	out_config[exp]['virus_name'] = virus
	out_config[exp]['virus_fasta'] = in_config[exp]['viruses'][virus]
	
	with open(out, 'w') as outfile:
		yaml.dump(out_config, outfile)


if __name__ == "__main__":
	main(sys.argv)
