#!/usr/bin/env python3

# usage: python3 make_isling_config.py <sim_config> <experiment_name> <sample_name> <output_config>


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
#  host_prefix: "test/references/test_human.fa"
#  virus_name: "rep68"
#  virus_prefix: "test/references/test_AAV.fa"
#  dedup: False
#  clip-cutoff: 20
#  min-mapq: 10
#  cigar-tol: 3
#  post:
#    - filter
#    - dedup
#  merge-dist: 1000
#  min-n-merge: 1
#  split: 5

# assume that there was only one host, virus and mean fragment length specified in the simulation config file

import sys
import yaml
import pdb
import os
import pandas as pd

def main(args):

	# get inputs
	config, exp, name, out = args[1:]
	
	# import simulation config yaml
	with open(config, 'r') as stream:
		try:
			in_config = yaml.safe_load(stream)
		except yaml.YAMLError as exc:
			print(exc)	

	# check that experiment specified is acutally in simulation config
	assert exp in in_config
	
	# apply global options to in config
	if 'global' in in_config:		
		# get default (global) options
		default = in_config.pop('global')
		for dataset in in_config:
			for key in default:
				if key not in in_config[dataset]:
					in_config[dataset][key] = default[key]
					
	# import simulation dataframe from results dir
	df_dir = os.path.join(in_config['coverage']['out_directory'], exp, "simulation_summary.tsv")
	in_df = pd.read_csv(df_dir, delimiter='\t')
	row_idx = in_df.index[in_df['unique'] == f"{exp}__{name}"][0]
	
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
	out_config[exp]['bam_suffix'] = '.sorted.bam'
	out_config[exp]['R1_suffix'] = '1.fq'
	out_config[exp]['R2_suffix'] = '2.fq'
	out_config[exp]['read1-adapt'] = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
	out_config[exp]['read2-adapt'] = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
	out_config[exp]['bwa-mem'] = "-A 1 -B 2 -O 6,6 -E 1,1 -L 0,0 -T 10 -h 200"
	
	# get read folder
	out_config[exp]['read_folder'] = os.path.join(in_config[exp]['out_directory'], exp, 'sim_reads')
	
	# set output directory
	out_config[exp]['out_dir'] = os.path.join(in_config[exp]['out_directory'], exp, "isling", name)	
	
	# set sample
	out_config[exp]['samples'] = [name]
	
	# host
	host = str(in_df.loc[row_idx, 'host_name'])
	out_config[exp]['host_name'] = host
	out_config[exp]['host_prefix'] = os.path.join(in_config[exp]['out_directory'], "references", host, host)
	out_config[exp]['host_prefix'] = os.path.normpath(out_config[exp]['host_prefix'])
	
	# virus
	virus = str(in_df.loc[row_idx, 'virus_name'])
	out_config[exp]['virus_name'] = virus
	out_config[exp]['virus_prefix'] = os.path.join(in_config[exp]['out_directory'], "references", virus, virus)
	out_config[exp]['virus_prefix'] = os.path.normpath(out_config[exp]['virus_prefix'])
	
	# mean fragment length
	out_config[exp]['mean-frag-len'] = float(in_df.loc[row_idx,'frag_len'])
	
	# cpus for alignment
	out_config[exp]['align-cpus'] = 20
	
	# number of parts for splitting - make this approximatley the fold coverage
	# but no more than 20 and no less than 1
	out_config[exp]['split'] = int(in_df.loc[row_idx,'fcov'])
	out_config[exp]['split'] = min(out_config[exp]['split'], 20)
	out_config[exp]['split'] = max(out_config[exp]['split'], 1)
	
	
	with open(out, 'w') as outfile:
		yaml.dump(out_config, outfile)
		
	print(f"saved output config to {outfile}")


if __name__ == "__main__":
	main(sys.argv)
