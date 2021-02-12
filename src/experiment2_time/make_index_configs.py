#!/usr/bin/env python3

# usage: python3 make_index_configs.py <sim_config> <vifi_data_repo_path> <isling_config> <other_tools_config> 


# given simulation config file, make config file for 
# indexing isling and other tools references on one sample from one experiment

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

# config file for other tools loooks like this:
#test:
#   read_directory: "data/reads"
#   out_directory: "test/out/"
#   R1_suffix: "1.fq"
#   R2_suffix: "2.fq"
#   read1-adapt: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
#   read2-adapt: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
#   analysis_host: 
#       "human" : "data/references/test_human.fa"
#   analysis_virus:
#       "AAV" :  "data/references/test_AAV.fa"
#   polyidus_params:
#     trim: True
#     aligner: "bowtie2"
#   vifi_params:
#     trim: True
#     host_info:
#       human: # fasta will be used from general analysis_hosts
#         mappability: "../vifi-test/data_repo/GRCh38/hg38full_k35_noMM.mappability.bedgraph"
#         mappability_exclude: "../data/references/GRCh38/ENCFF356LFX.bed"
#         genes: "../data/references/GRCh38/hg38.gencode.v35.annotation.genes.gff3"
#         exons: "../data/references/GRCh38/hg38.gencode.v35.annotation.exons.gff3"
#         oncogenes: "../data/references/GRCh38/hg38.gencode.v35.annotation.genes.gff3"
#         centromeres: "../data/references/GRCh38/centromeres.bed"
#         conserved_regions: "../data/references/GRCh38/exclude.cnvnator_100bp.GRCh38.20170403.bed"
#         segdup: "../data/references/GRCh38/genomicSuperDups.bed"  
#   seeksv_params:
#     trim: True
#     dedup: False
#   vseq_toolkit_params:
#     qua: 20
#     lenPer: 50
#     vecVecFusion: 'true'
#     stringencyVec: 'high'
#     UMthresholdVec: 0.95
#     minMapSpanVec: 20
#     distVecVec: 10
#     opVecVec: 5
#     idenVecVec: 95
#     stringencyVecGen: 'high'
#     UMthresholdVecGen: 0.95
#     minMapSpanVecGen: 20
#     distVecGen: 10
#     opVecGen: 5
#     idenVecGen: 95
#     clusterRange: 3
#     host_table:
#       human: "$VSeqToolkit/testDir/testReferenceIndex/refSeqUCSCTablehg38.txt"


# assume that all experiments have the same host and virus

import sys
import yaml
import pdb
import os

def main(args):

	# get inputs
	in_config_path, vifi_repo, isling_config_path, tools_config_path = args[1:]
	
	# import simulation config yaml
	with open(in_config_path, 'r') as stream:
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
	
	make_isling_config(in_config, isling_config_path)
	make_tools_config(in_config, tools_config_path, vifi_repo)
	
def make_tools_config(in_config, tools_config_path, vifi_repo):
	exp = list(in_config.keys())[0]
	host = list(in_config[exp]['hosts'].keys())[0]
	virus = list(in_config[exp]['viruses'].keys())[0]
	
	seeksv_params = {
		'trim': False,
		'dedup': False
	}
	polyidus_params = {
		'trim': False,
		'aligner': 'bowtie2'
	}
	vifi_params = get_vifi_params(vifi_repo, host)
	vseq_params = {
		'qua': 20,
		'lenPer' : 50,
		'mode' : 'default',
		'vecVecFusion': 'true',
		'stringencyVec': 'medium',
		'UMthresholdVec' : 0.95,
		'minMapSpanVec': 20,
		'distVecVec' : 10,
		'opVecVec' : 5,
		'idenVecVec' : 95,
		'stringencyVecGen' : 'medium',
		'UMthresholdVecGen' : 0.95,
		'minMapSpanVecGen' : 20,
		'distVecGen' : 10,
		'opVecGen': 5,
		'idenVecGen' : 95,
		'clusterRange': 3,
		'host_table' : {host : '$VSeqToolkit/testDir/testReferenceIndex/refSeqUCSCTablehg38.txt'}
	}
	
	out_config = {
		exp: {
			'read_directory': os.path.join(in_config[exp]['out_directory'], exp, 'sim_reads'),
			'out_directory': os.path.normpath(in_config[exp]['out_directory']),
			'R1_suffix': '1.fq',
			'R2_suffix': '2.fq',
			'read1-adapt' : "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
			'read2-adapt' : "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",	
			'analysis_host' : {host: in_config[exp]['hosts'][host]},
			'analysis_virus' : {virus: in_config[exp]['viruses'][virus]},
			'polyidus_params': polyidus_params,
			'vifi_params': vifi_params,
			'seeksv_params': seeksv_params,
			'vseq_toolkit_params': vseq_params
		}
	}
	
	with open(tools_config_path, 'w') as outfile:
		yaml.dump(out_config, outfile)

	print(f"saved other tools config to {tools_config_path}")
	
def get_vifi_params(vifi_repo_path, host_name):
	
	#get references in repo
	refdir = os.path.abspath(os.path.join(vifi_repo_path,"GRCh38"))
	ref_files = {}
	with open(os.path.join(refdir, "file_list.txt"), 'r') as handle:
		for line in handle:
			f = line.strip().split()
			ref_files[f[0]] = f[1]
			
	return {
		'trim': False,
		'host_info': {
			host_name : {
				'mappability': os.path.join(refdir, ref_files['duke35_filename']),
				'mappability_exclude': os.path.join(refdir, ref_files['mapability_exclude_filename']),
				'genes': os.path.join(refdir, ref_files['gene_filename']),
				'exons': os.path.join(refdir, ref_files['gene_filename']),
				'oncogenes': os.path.join(refdir, ref_files['oncogene_filename']),
				'centromeres': os.path.join(refdir, ref_files['centromere_filename']),
				'conserved_regions': os.path.join(refdir, ref_files['conserved_regions_filename']),
				'segdup': os.path.join(refdir, ref_files['segdup_filename']),
			}
		}
	}

def make_isling_config(in_config, isling_config_path):
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
	out_config[exp]['out_dir'] = os.path.normpath(in_config[exp]['out_directory'])	
	
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
	
	with open(isling_config_path, 'w') as outfile:
		yaml.dump(out_config, outfile)

	print(f"saved isling config to {isling_config_path}")

if __name__ == "__main__":
	main(sys.argv)
