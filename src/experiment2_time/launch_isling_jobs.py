#!/usr/bin/env python3

# run isling on each pair of fastq files from simulation
# collect runtime from each and save to file

# usage: python3 launch_isling_jobs.py <config file> <container_name> <config_script>

# collect information from time 
#           %Uuser %Ssystem %Eelapsed %PCPU (%Xtext+%Ddata %Mmax)k
#           %Iinputs+%Ooutputs (%Fmajor+%Rminor)pagefaults %Wswaps
#
#       %U     Total number of CPU-seconds that the process spent in user mode.
#       %S     Total number of CPU-seconds that the process spent in kernel mode.
#       %E     Elapsed real time (in [hours:]minutes:seconds).
#       %P     Percentage of the CPU that this job got, computed as (%U + %S) / %E.
#       %X     Average size of the process's shared text space, in Kbytes.
#       %D     Average size of the process's unshared data area, in Kbytes.
#       %M     Maximum resident set size of the process during its lifetime, in Kbytes.
#       %I     Number of filesystem inputs by the process.
#       %O     Number of filesystem outputs by the process.
#       %F     Number of major page faults that occurred while the process was running.  These are faults where the page
#              has to be read in from disk.
#       %R     Number  of  minor,  or recoverable, page faults.  These are faults for pages that are not valid but which
#              have not yet been claimed by other virtual pages.  Thus the data in the page is still valid but the  sys-
#              tem tables must be updated.
#       %W     Number of times the process was swapped out of main memory.


import sys
import os
import yaml
import pdb
import glob
import csv
import subprocess

replicates = 3

def main(argv):
	
	# get args
	config_path, container, config_script = argv[1:]
	
	# import config file
	sim_config = import_yaml(config_path)
	
	
	results = []
	# get sample in each dir in config file
	for exp in sim_config:
		samples = get_samples(exp, sim_config[exp])
		
		# run isling for each sample
		for samp in samples:
			results += run_isling(sim_config, exp, samp, config_path, config_script, container, replicates)
			
		# write output
		write_output(sim_config, exp, results)	

	
def write_output(sim_config, exp, results):
	# write results
	outfile = os.path.join(sim_config[exp]['out_directory'], exp, "isling", f"{exp}_runtime.tsv")
	with open(outfile, 'w', newline='') as outhandle:
		fieldnames = ['tool', 'dataset', 'sample', 'replicate', 
									'user_time', 'system_time', 'elapsed_time', 'CPU', 'shared_text', 'unshared_data',
									'max_rss', 'fs_inputs', 'fs_outputs', 'major_page_faults', 'minor_page_faults', 'swaps'
		]
		writer = csv.DictWriter(outhandle, fieldnames=fieldnames, delimiter='\t')
			
		writer.writeheader()
		for r in results:
			writer.writerow(r)
			
			
def run_isling(sim_config, exp, sample, sim_config_path, config_script, container, reps):

	isling_config = os.path.join(sim_config[exp]['out_directory'], exp, 
																'isling', sample, f"{exp}_{sample}.yml")	
	
	# make directory for config file
	isling_dir = os.path.dirname(isling_config)
	a = subprocess.run(['mkdir', '-p', isling_dir])
	a.check_returncode()
	
	# make config file
	args = ['python3', config_script, sim_config_path, exp, sample, isling_config]
	a = subprocess.run(args)
	a.check_returncode()
	
	# run isling
	config = import_yaml(isling_config)
	# host prefix dir
	host_prefix = os.path.realpath(os.path.dirname(config[exp]['host_prefix']))
	# virus prefix dir
	virus_prefix = os.path.realpath(os.path.dirname(config[exp]['virus_prefix']))
	# read dir
	read_dir = os.path.realpath(config[exp]['read_folder'])
	# output dir
	out_dir = os.path.realpath(config[exp]['out_dir'])
	
	srun = ['srun', '--exclusive', '-c20', '-n1', '--mem', '128gb']
	time = ['/usr/bin/time']
	sing = ['singularity', 'exec',
					'-B', host_prefix,
					'-B', virus_prefix,
					'-B', read_dir,
					'-B', out_dir,
					container
				]
	isling = ['snakemake', '--jobs', '100', '--configfile', isling_config, '--forceall']
		
	args =  srun + time + sing + isling
	results = []
	for i in range(reps):
		a = subprocess.run(args, capture_output=True, text=True)
		a.check_returncode()
	
		res = collect_output(a.stderr)
		res['tool'] = 'isling'
		res['dataset'] = exp
		res['sample'] = sample
		res['replicate'] = i
		results.append(res)
	
	return results



def collect_output(stdout):
#   %Uuser %Ssystem %Eelapsed %PCPU (%Xtext+%Ddata %Mmax)k
#   %Iinputs+%Ooutputs (%Fmajor+%Rminor)pagefaults %Wswaps
		stdout = stdout.split('\n')
		line1 = stdout[-3].split()
		line2 = stdout[-2].split()
		results = {
			'user_time' : line1[0][:-4],
			'system_time': line1[1][:-6],
			'elapsed_time': line1[2][:-7],
			'CPU': line1[3][:-3],
			'shared_text': line1[4].split('+')[0][1:-7],
			'unshared_data': line1[4].split('+')[1][:-7],			
			'max_rss': line1[5][:-13],		
			'fs_inputs': line2[0].split('+')[0][:-6],
			'fs_outputs': line2[0].split('+')[1][:-7],			
			'major_page_faults': line2[1].split('+')[0][1:-5],			
			'minor_page_faults': line2[1].split('+')[1][:-16],			
			'swaps': line2[2][:-5],		
		}
		
		return results
		
	
def get_samples(exp, config_exp):

	readdir = os.path.join(config_exp['out_directory'], exp, 'sim_reads')	
	suffix = "*1.fq"
	
	samples = glob.glob(os.path.join(readdir, suffix))
	samples = [s[len(readdir)+1:-(len(suffix)-1)] for s in samples]
	return samples

def import_yaml(config):
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
	
	return in_config

if __name__ == "__main__":
	main(sys.argv)
