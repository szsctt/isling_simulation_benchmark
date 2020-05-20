#!/bin/bash
#
#SBATCH --job-name=integrate_test		#name of job
#SBATCH --mail-type=BEGIN,END,FAIL		#email the events (NONE,BEGIN,END,FAIL,ALL)
#SBATCH --output=../../slogs/chr20/overlapwhole/300int_rep4MATCH.txt #standard error and output log
#SBATCH --mail-user=Susie.Grigson@cisro.au	#where to send mail
#SBATCH --nodes=1				#number of nodes 
#SBATCH --ntasks=1				#number of tasks 
#SBATCH --cpus-per-task=1			#number of CPUs per task
#SBATCH	--mem=18GB 				#memory per NODE 
#SBATCH	--time=01:30:00				#time limit hrs:min:sec
#SBATCH	--error=../../slogs/chr20/overlapwhole/300int_rep4MATCH_error.txt 		#error destination 

pwd; hostname; date

echo "Running on $SLURM_JOB_NUM_NODES nodes with $SLURM_NTASKS tasks, each with $SLURM_CPUS_PER_TASK cores."

echo "Batch is $SBATCH_BATCH, job is $SBATCH_JOB_NAME"

srun python3 read_integrations.py --sam ../../script_outputs/chr20/overlapwhole/300int/rep_4/paired_dat.sam --host_ints ../../script_outputs/chr20/overlapwhole/300int/rep_4/host_integrations.csv --save ../../script_outputs/chr20/overlapwhole/300int/rep_4/all_reads.csv --viral_only ../../script_outputs/chr20/overlapwhole/300int/rep_4/viral_only.csv 

wait 

date 
