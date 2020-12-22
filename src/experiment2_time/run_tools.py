#!/usr/bin/env python3

# run isling on each pair of fastq files from simulation
# collect runtime from each and save to file

# usage: python3 run_tools.py <config file> <isling_sif> <isling_config_script>
# <seeksv_sif> <seeksv_script> <polyidus_sif>

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
import multiprocessing as mp
import multiprocessing.pool as mp_pool
import functools
import time
import random

replicates = 3
seeksv_threads=20
srun_args = ['srun', '--exclusive', '-c20', '-n1', '--mem', '128gb', '--time', '24:00:00']
srun_debug = ['srun']
#srun_args = srun_debug
time_args = ['/usr/bin/time']
fieldnames = ['tool', 'dataset', 'sample', 'replicate', 'exit_value',
							'user_time', 'system_time', 'elapsed_time', 'CPU', 'shared_text', 'unshared_data',
							'max_rss', 'fs_inputs', 'fs_outputs', 'major_page_faults', 'minor_page_faults', 'swaps'
		]

def main(argv):
	
	# get args
	[config_path, isling_container, config_script,
	seeksv_container, seeksv_script, polyidus_container,
	vifi_data_repo, vifi_container] = argv[1:]
	
	# import config file
	sim_config = import_yaml(config_path)

	# get sample in each dir in config file
	run = functools.partial(run_experiment, sim_config=sim_config, 
														config_path=config_path, config_script=config_script, 
														seeksv_script=seeksv_script, isling_container=isling_container,
														seeksv_container=seeksv_container, polyidus_container=polyidus_container,
														vifi_data_repo=vifi_data_repo, vifi_container=vifi_container,
														reps=replicates
													)
	procs = []
	with mp_pool.ThreadPool(2) as pool:
		for exp in sim_config:
		
			# start process
			procs.append(pool.apply_async(run, (exp, )))
			#result = run(exp)

		
		[p.get() for p in procs]

def run_experiment(exp, sim_config, config_path, config_script, seeksv_script, isling_container, seeksv_container, polyidus_container, vifi_data_repo, vifi_container, reps):

	# write header for experiment
	outfile = os.path.join(sim_config[exp]['out_directory'], exp, f"{exp}_runtime.tsv")
	if not os.path.isfile(outfile):
		with open(outfile, 'w', newline='') as outhandle:
			writer = csv.DictWriter(outhandle, fieldnames=fieldnames, delimiter='\t')
			writer.writeheader()
	
	# make partial functions
	m = mp.Manager()
	lock = m.Lock()
	run_isling_partial = functools.partial(run_isling, sim_config=sim_config, 
																					sim_config_path=config_path, config_script=config_script,		
																				 container=isling_container, reps=reps, outfile=outfile, 
																				 lock=lock)
	run_seeksv_partial = functools.partial(run_seeksv, sim_config=sim_config, seeksv_script=seeksv_script, 
																					container=seeksv_container, reps=reps, outfile=outfile, \
																					lock=lock)
	run_polyidus_partial = functools.partial(run_polyidus, sim_config=sim_config, 
																				 	container=polyidus_container, reps=reps, 
																				 	outfile=outfile, lock=lock)	
	run_vifi_partial = functools.partial(run_vifi, sim_config=sim_config,vifi_data_repo = vifi_data_repo,
																				 	container=vifi_container, reps=reps, 
																				 	outfile=outfile, lock=lock, threads=seeksv_threads)	
	# get samples for this experiment
	samples = get_samples(exp, sim_config[exp])

	with mp.Pool(10) as pool:

		procs = []
		for samp in samples:
			#procs.append(pool.apply_async(run_isling_partial, (exp, samp)))
			#procs.append(pool.apply_async(run_seeksv_partial, (exp, samp)))
			#procs.append(pool.apply_async(run_polyidus_partial, (exp, samp)))
					
			run_vifi_partial(exp, samp)
		# get all results
		[p.get() for p in procs]
		

def run_polyidus(exp, sample, sim_config, container, reps, outfile, lock):
	print(f"running polyidus on experiment {exp}, sample {sample} at time {time.strftime('%a, %d %b %Y %H:%M:%S', time.localtime())}")

	# make directory for output
	poly_dir = os.path.join(sim_config[exp]['out_directory'], exp, 
																'polyidus', sample)
	a = subprocess.run(['mkdir', '-p', poly_dir])
	a.check_returncode()	
	
	# collect inputs
	r1 = os.path.join(sim_config[exp]['out_directory'], exp, 'sim_reads', f"{sample}1.fq")
	r2 = os.path.join(sim_config[exp]['out_directory'], exp, 'sim_reads', f"{sample}2.fq")
	virus = list(sim_config[exp]['viruses'].keys())[0]
	virus_prefix = os.path.join(sim_config[exp]['out_directory'], 'references', virus, virus)
	host = list(sim_config[exp]['hosts'].keys())[0]
	host_prefix = os.path.join(sim_config[exp]['out_directory'], 'references', host, host)
	
	# arguments
	sing_args = ['singularity', 'exec', 
								'-B', os.path.abspath(os.path.dirname(r1)),
								'-B', os.path.abspath(os.path.dirname(host_prefix)),
								'-B', os.path.abspath(os.path.dirname(virus_prefix)),
								'-B', os.path.abspath(poly_dir),
								container]		
	poly_args = ['python3', '/usr/src/app/src/polyidus.py', host_prefix, virus_prefix,
								'--fastq', r1, r2, '--outdir', poly_dir]


	# this one can't be parallel because they have the same output files
	for i in range(reps):
		if check_already_run(lock, outfile, 'seeksv', exp, sample, i):
			continue
		args =  srun_args + ['--job-name', f'seeksv.{exp}.{sample}.{i}'] + time_args + sing_args + poly_args
		polyidus = subprocess.run(args, capture_output=True, text=True)
		collect_output(polyidus, 'polyidus', exp, sample, i, outfile, lock)


def run_vifi(exp, sample, sim_config, vifi_data_repo, container, reps, outfile, lock, threads):
	print(f"running vifi on experiment {exp}, sample {sample} at time {time.strftime('%a, %d %b %Y %H:%M:%S', time.localtime())}")

	# make directory for output
	vifi_dir = os.path.join(sim_config[exp]['out_directory'], exp, 
																'vifi', sample)
	a = subprocess.run(['mkdir', '-p', vifi_dir])
	a.check_returncode()
	
	os.environ['AA_DATA_REPO'] = vifi_data_repo
	os.environ['REFERENCE_REPO'] = os.path.join(sim_config[exp]['out_directory'], "vifi_refs", "data_repo")	
	os.environ['VIFI_DIR'] ='/home/ViFi'
	
	# make chromosome list
	host = list(sim_config[exp]['hosts'].keys())[0]
	host_fasta = sim_config[exp]['hosts'][host]

	host_fai = host_fasta + ".fai"
	chr_list = os.path.join(os.path.dirname(host_fasta), f"{host}_chrlist.txt")
	# if the chromosome list doesn't exist, make it
	with lock:
		if not os.path.isfile(host_fai):
			args = ['singularity', 'exec', '-B', os.path.abspath(os.path.dirname(host_fai)), 
								container, 'samtools', 'faidx', host_fasta]
			a = subprocess.run(args)
			a.check_returncode()
		if not os.path.isfile(chr_list):
			awk = 'BEGIN {{ ORS = " " }} {{a[$1]}} END {{for (i in a) print i}}'
			args = ['awk', awk, host_fai]
			a = subprocess.run(args, capture_output=True, text=True)
			a.check_returncode()
			with open(chr_list, 'w') as handle:
				handle.write(a.stdout)
	
	# collect other arguments
	virus = list(sim_config[exp]['viruses'].keys())[0]
	r1 = os.path.join(sim_config[exp]['out_directory'], exp, 'sim_reads', f"{sample}1.fq")
	r2 = os.path.join(sim_config[exp]['out_directory'], exp, 'sim_reads', f"{sample}2.fq")
	ref = os.path.join(sim_config[exp]['out_directory'], 'vifi_refs', 'data', virus, f"{host}_{virus}.fas")
		
	sing_args = ['singularity', 'exec', '-B', os.path.abspath(os.path.dirname(r1)),
								'-B', os.path.abspath(os.path.dirname(ref)),
								'-B', os.path.abspath(os.path.dirname(chr_list)),
								'-B', os.path.abspath(vifi_data_repo),
								'-B', os.path.abspath(os.path.join(sim_config[exp]['out_directory'], 
																			"vifi_refs", "data")),
							container]	
	
	vifi_args = ['/usr/bin/python', '/home/ViFi/scripts/run_vifi.py', '-f', r1, '-r', r2,
								 '-c', str(threads), '--reference', ref,
								'-v', virus, '-o', vifi_dir, '-d', 'True', '-C', chr_list]
		
	for i in range(reps):
		if check_already_run(lock, outfile, 'vifi', exp, sample, i):
			continue
		args =  srun_args + ['--job-name', f'vifi.{exp}.{sample}.{i}'] + time_args + sing_args + vifi_args
		seeksv = subprocess.run(args, capture_output=True, text=True)
		collect_output(seeksv, 'vifi', exp, sample, i, outfile, lock)	

def run_seeksv(exp, sample, sim_config, seeksv_script, container, reps, outfile, lock):
	print(f"running seeksv on experiment {exp}, sample {sample} at time {time.strftime('%a, %d %b %Y %H:%M:%S', time.localtime())}")
	
	# make directory for output
	seeksv_dir = os.path.join(sim_config[exp]['out_directory'], exp, 
																'seeksv', sample)
	a = subprocess.run(['mkdir', '-p', seeksv_dir])
	a.check_returncode()
																	
	# collect inputs for seeksv script
	r1 = os.path.join(sim_config[exp]['out_directory'], exp, 'sim_reads', f"{sample}1.fq")
	r2 = os.path.join(sim_config[exp]['out_directory'], exp, 'sim_reads', f"{sample}2.fq")
	virus = list(sim_config[exp]['viruses'].keys())[0]
	host = list(sim_config[exp]['hosts'].keys())[0]
	prefix = os.path.join(sim_config[exp]['out_directory'], 'seeksv_refs', virus, f"{host}_{virus}.fas")
	threads = seeksv_threads
	
	sing_args = ['singularity', 'exec', 
								'-B', os.path.abspath(os.path.dirname(r1)),
								'-B', os.path.abspath(os.path.dirname(prefix)),
								'-B', os.path.abspath(seeksv_dir),
								container]
	seeksv_args = ['bash', seeksv_script, r1, r2, prefix, str(threads), seeksv_dir]

	# this one can't be parallel because they have the same output files
	for i in range(reps):
		if check_already_run(lock, outfile, 'seeksv', exp, sample, i):
			continue
		args =  srun_args + ['--job-name', f'seeksv.{exp}.{sample}.{i}'] + time_args + sing_args + seeksv_args
		seeksv = subprocess.run(args, capture_output=True, text=True)
		collect_output(seeksv, 'seeksv', exp, sample, i, outfile, lock)

			
def run_isling(exp, sample, sim_config, sim_config_path, config_script, container, reps, outfile, lock):

	print(f"running isling on experiment {exp}, sample {sample} at time {time.strftime('%a, %d %b %Y %H:%M:%S ', time.localtime())}")
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
	

	sing_args = ['singularity', 'exec',
					'-B', host_prefix,
					'-B', virus_prefix,
					'-B', read_dir,
					'-B', out_dir,
					container
				]
	isling_args = ['snakemake', '--jobs', '100', '--configfile', isling_config, '--forceall']
		
	# this one can't be parallel because they have the same output files
	# and snakemake will chuck a hissy 
	for i in range(reps):
		if check_already_run(lock, outfile, 'seeksv', exp, sample, i):
			continue
		args =  srun_args + ['--job-name', f'isling.{exp}.{sample}.{i}'] + time_args + sing_args + isling_args
		isling = subprocess.run(args, capture_output=True, text=True)
		collect_output(isling, 'isling', exp, sample, i, outfile, lock)


def collect_output(proc, tool, dataset, sample, rep, outfile, lock):
#   %Uuser %Ssystem %Eelapsed %PCPU (%Xtext+%Ddata %Mmax)k
#   %Iinputs+%Ooutputs (%Fmajor+%Rminor)pagefaults %Wswaps
		
	if proc.returncode == 0:
		stdout = proc.stderr.split('\n')
		line1 = stdout[-3].split()
		line2 = stdout[-2].split()
	else:
		print(proc.stderr)
		line1 = ['user', 'system', 'elapsed', 'CPU', '(text+data', 'max)k']
		line2 = ['inputs+outputs', '(major+minor)pagefaults', 'swaps'] * 10

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
			'tool' : tool,
			'dataset': dataset,
			'sample': sample,
			'replicate': rep,
			'exit_value': proc.returncode
					
	}
	
	with lock, open(outfile, 'a', newline='') as outhandle:	
		writer = csv.DictWriter(outhandle, fieldnames=fieldnames, delimiter='\t')
		writer.writerow(results)
			

def check_already_run(lock, outfile, tool, experiment, sample, rep):
	
	with lock, open(outfile, 'r', newline='') as outhandle:
		reader = csv.DictReader(outhandle, fieldnames=fieldnames, delimiter='\t')
		for row in reader:
			if row['tool'] != tool:
				continue
			if row['dataset'] != experiment:
				continue
			if row['sample'] != sample:
				continue
			if row['replicate'] != str(rep):
				continue
			return True
			
	return False
	

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
