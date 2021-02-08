#!/usr/bin/env python3

# run isling on each pair of fastq files from simulation
# collect runtime from each and save to file

# usage: python3 run_tools.py <config file> <isling_sif> <isling_config_script>
# <seeksv_sif> <seeksv_script> <polyidus_sif> <vseq_toolkit_sif>

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
import argparse
max_time = 86400 # max time in seconds for each run of each tool
replicates = 3
retries = 1
seeksv_threads = 20
srun_args = ['srun', '--exclusive', '-c20', '-n1', '--mem', '128gb', '--time', '24:00:00']
srun_debug = ['srun']
#srun_args = srun_debug
time_args = ['/usr/bin/time']
fieldnames = ['tool', 'dataset', 'sample', 'replicate', 'exit_value',
							'user_time', 'system_time', 'elapsed_time', 'CPU', 'shared_text', 'unshared_data',
							'max_rss', 'fs_inputs', 'fs_outputs', 'major_page_faults', 'minor_page_faults', 'swaps', 'command'
		]

def main(argv):

	# get args
	parser = argparse.ArgumentParser(description='run tools on simulated data and collect runtimes')
	parser.add_argument('--sim-config', '-c', help='config from simulation')
	parser.add_argument('--isling-sif', help='path to isling sif file')
	parser.add_argument('--isling-config-script', help='path to script that is used to create config files for isiling')
	parser.add_argument('--seeksv-sif', help='path to seeksv sif file')
	parser.add_argument('--seeksv-script', help='bash script for running seeksv on one sample')
	parser.add_argument('--polyidus-sif', help='path to polyidus sif file')
	parser.add_argument('--vifi-sif', help='path to vifi sif file')	
	parser.add_argument('--vifi-data-repo', help='path to vifi data repo')	
	parser.add_argument('--vseq-toolkit-sif', help='path to vifi sif file')
	parser.add_argument('--parallel', action='store_true', help='run jobs in parallel on the cluster?')
	parser.add_argument('--replicates', help='number of replicates to perform', default=replicates, type=int)	
	parser.add_argument('--cores', help='number of cores to use', default=16, type=int)
	args = parser.parse_args(argv[1:])
	
	# import config file
	sim_config = import_yaml(args.sim_config)

	# get sample in each dir in config file
	run = functools.partial(run_experiment, sim_config=sim_config, 
														config_path=args.sim_config, config_script=args.isling_config_script, 
														seeksv_script=args.seeksv_script, isling_container=args.isling_sif,
														seeksv_container=args.seeksv_sif, polyidus_container=args.polyidus_sif,
														vifi_data_repo=args.vifi_data_repo, vifi_container=args.vifi_sif,
														reps=args.replicates, parallel=args.parallel, cores=args.cores
													)
	procs = []
	with mp_pool.ThreadPool(2) as pool:
		for exp in sim_config:
		
			# start processes
			if args.parallel:
				procs.append(pool.apply_async(run, (exp, )))
			else:
				result = run(exp)

		if args.parallel:
			[p.get() for p in procs]

def run_experiment(exp, sim_config, config_path, config_script, seeksv_script, isling_container, seeksv_container, polyidus_container, vifi_data_repo, vifi_container, reps, parallel, cores):

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
																				 lock=lock, retries=retries, parallel=parallel, threads=cores)
	run_seeksv_partial = functools.partial(run_seeksv, sim_config=sim_config, seeksv_script=seeksv_script, 
																					container=seeksv_container, reps=reps, outfile=outfile, \
																					lock=lock, retries=retries, parallel=parallel, threads=cores)
	run_polyidus_partial = functools.partial(run_polyidus, sim_config=sim_config, 
																				 	container=polyidus_container, reps=reps, 
																				 	outfile=outfile, lock=lock, retries=retries, parallel=parallel)	
	run_vifi_partial = functools.partial(run_vifi, sim_config=sim_config,vifi_data_repo = vifi_data_repo,
																				 	container=vifi_container, reps=reps, 
																				 	outfile=outfile, lock=lock, threads=cores,
																				 	retries=retries, parallel=parallel)	
	# get samples for this experiment
	samples = get_samples(exp, sim_config[exp])

	with mp.Pool(10) as pool:

		procs = []
		for samp in samples:
			if parallel:
				procs.append(pool.apply_async(run_isling_partial, (exp, samp)))
				#procs.append(pool.apply_async(run_seeksv_partial, (exp, samp)))
				#procs.append(pool.apply_async(run_polyidus_partial, (exp, samp)))
				#procs.append(pool.apply_async(run_vifi_partial, (exp, samp)))					
			else:
				run_isling_partial(exp, samp)
				run_seeksv_partial(exp, samp)
				run_polyidus_partial(exp, samp)
				run_vifi_partial(exp, samp)
				
		#get all results
		if parallel:
			[p.get() for p in procs]

def run_vseq_toolkit(exp, sample, sim_config, container, reps, outfile, lock, retries, parallel):

	# make directory for output
	vseq_dir = os.path.join(sim_config[exp]['out_directory'], exp, 
																'vseq_toolkit', sample)
	a = subprocess.run(['mkdir', '-p', vseq_dir])
	a.check_returncode()	
	
	# collect inputs
	r1 = os.path.join(sim_config[exp]['out_directory'], exp, 'sim_reads', f"{sample}1.fq")
	r2 = os.path.join(sim_config[exp]['out_directory'], exp, 'sim_reads', f"{sample}2.fq")
	virus = list(sim_config[exp]['viruses'].keys())[0]
	virus_prefix = os.path.join(sim_config[exp]['out_directory'], 'references', virus, virus)
	host = list(sim_config[exp]['hosts'].keys())[0]
	host_prefix = os.path.join(sim_config[exp]['out_directory'], 'references', host, host)
	
	# write config file
	config = os.path.join(vseq_dir, "config.txt")
	with open(config, 'w') as handle:
		handle.write(f"file1= {r1}\n")
		handle.write(f"file2= {r2}\n")
		handle.write(f"outDir= {vseq_dir}\n")
		handle.write("bin= $VSeqToolkit/scripts/\n")
		handle.write("qua=20\n")
		handle.write("lenPer=50\n")
		handle.write("adapter1=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA\n")
		handle.write("adapter2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT\n")
		handle.write("trimmer= $VSeqToolkit/thirdPartyTools/skewer\n")
		handle.write("aligner= $VSeqToolkit/thirdPartyTools/bwa\n")
		handle.write("samtools= $VSeqToolkit/thirdPartyTools/samtools\n")
		handle.write("mode=default\n")
	
	
	# arguments
	sing_args = ['singularity', 'exec', 
								'-B', os.path.abspath(os.path.dirname(r1)),
								'-B', os.path.abspath(os.path.dirname(host_prefix)),
								'-B', os.path.abspath(os.path.dirname(virus_prefix)),
								'-B', os.path.abspath(poly_dir),
								container]		
	poly_args = ['perl', '/usr/src/app/src/polyidus.py', '-c', config]


	# run polyidus
	run_tool(lock, outfile, poly_dir, 'polyidus', exp, sample, reps, retries, parallel, sing_args, poly_args)


def run_polyidus(exp, sample, sim_config, container, reps, outfile, lock, retries, parallel):

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


	# run polyidus
	run_tool(lock, outfile, poly_dir, 'polyidus', exp, sample, reps, retries, parallel, sing_args, poly_args)

def run_vifi(exp, sample, sim_config, vifi_data_repo, container, reps, outfile, lock, threads, retries, parallel):

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

	run_tool(lock, outfile, vifi_dir, 'vifi', exp, sample, reps, retries, parallel, sing_args, vifi_args)

def run_seeksv(exp, sample, sim_config, seeksv_script, container, reps, outfile, lock, retries, parallel, threads):
	
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
	
	sing_args = ['singularity', 'exec', 
								'-B', os.path.abspath(os.path.dirname(r1)),
								'-B', os.path.abspath(os.path.dirname(prefix)),
								'-B', os.path.abspath(seeksv_dir),
								container]
	seeksv_args = ['bash', seeksv_script, r1, r2, prefix, str(threads), seeksv_dir]

	run_tool(lock, outfile, seeksv_dir, 'seeksv', exp, sample, reps, retries, parallel, sing_args, seeksv_args)
			
def run_isling(exp, sample, sim_config, sim_config_path, config_script, container, reps, outfile, lock, retries, parallel, threads):

	isling_config = os.path.join(sim_config[exp]['out_directory'], exp, 
																'isling', sample, f"{exp}_{sample}.yml")	
	
	# make directory for config file
	isling_dir = os.path.dirname(isling_config)
	a = subprocess.run(['mkdir', '-p', isling_dir])
	a.check_returncode()
	
	# make config file
	config_args = ['python3', config_script, sim_config_path, exp, sample, isling_config, str(threads)]
	a = subprocess.run(config_args)
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
		
	# run tool 
	run_tool(lock, outfile, isling_dir, 'isling', exp, sample, reps, retries, parallel, sing_args, isling_args)


def run_tool(lock, outfile, out_dir, tool, exp, sample, reps, retries, parallel, sing_args, tool_args):

	# this one can't be parallel because all runs with the same data have the same output files
	for i in range(reps):
		if check_already_run(lock, outfile, tool, exp, sample, i):
			print(f"skipping {tool} on experiment {exp}, sample {sample} at time {time.strftime('%a, %d %b %Y %H:%M:%S', time.localtime())}")
			continue
		run_args = time_args + sing_args + tool_args
		if parallel:
			run_args = srun_args + ['--job-name', f'{tool}.{exp}.{sample}.{i}'] + run_args
		for j in range(retries):
			log = os.path.join(out_dir, f"{sample}.run_rep{i}.try{j}.log")
			try:
				print(f"running try {j}, replicate {i} for {tool} on experiment {exp}, sample {sample} at time {time.strftime('%a, %d %b %Y %H:%M:%S', time.localtime())}")
				if parallel:
					tool_run = subprocess.run(run_args, capture_output=True, text=True)
				else:
					tool_run = subprocess.run(run_args, capture_output=True, text=True, timeout=max_time)  
				
			except subprocess.TimeoutExpired:
				print('timed out!')
				handle_error(tool_run, log, exp, sample, tool, i, j)
				break
			if tool_run.returncode == 0:
				write_log(tool_run, log)
				break
			handle_error(tool_run, log, exp, sample, tool, i, j)

		try:
			collect_output(tool_run, tool, exp, sample, i, outfile, lock, run_args)
		except NameError:
			write_blank_line('timeout', tool, exp, sample, i, outfile, lock, run_args)


def handle_error(tool_run, log_file, exp, sample, tool, rep, retry):
	print(f'{tool} failed while trying to process sample {sample} from experiment {exp} on try {retry} of replicate {rep}')
	print(f'see log file {os.path.abspath(log_file)} for more information\n')
	write_log(tool_run, log_file)

def write_log(tool_run, log_file):	
	output = f"stdout:\n{tool_run.stdout}\n\nstderr:\n{tool_run.stderr}\n"
	with open (log_file, 'w') as handle:
		handle.write(output)	

def write_blank_line(return_code, tool, dataset, sample, rep, outfile, lock, args):
	results = { 
		'user_time' : '',
		'system_time': '',
		'elapsed_time': '', 
		'CPU': '',  
		'shared_text': '',
		'unshared_data': '',
		'max_rss': '',
		'fs_inputs': '',
		'fs_outputs': '',
		'major_page_faults': '',
		'minor_page_faults': '',
		'swaps': '',
		'tool' : tool,
		'dataset': dataset,
		'sample': sample,
		'replicate': rep,  
		'exit_value': return_code, 
		'command': " ".join(args)
	}

	with lock, open(outfile, 'a', newline='') as outhandle:
		writer = csv.DictWriter(outhandle, fieldnames=fieldnames, delimiter='\t')
		writer.writerow(results)

def collect_output(proc, tool, dataset, sample, rep, outfile, lock, args):
#   %Uuser %Ssystem %Eelapsed %PCPU (%Xtext+%Ddata %Mmax)k
#   %Iinputs+%Ooutputs (%Fmajor+%Rminor)pagefaults %Wswaps
		
	if proc.returncode == 0:
		stdout = proc.stderr.split('\n')
		line1 = stdout[-3].split()
		line2 = stdout[-2].split()
	else:
		#print(proc.stderr)
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
			'exit_value': proc.returncode,
			'command': " ".join(args)
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
