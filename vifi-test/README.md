## Testing vifi

Tested ViFi on a few different datasets and with a few different parameters.  This tool seems mainly to be geared towards finding wild-type viruses that may or may not be like the viral references provided.  It uses a hidden markov model to decide if reads not aligned to the virus are similar enough that they might actually be viral.  This requires building such models - they are provided for hpv and hbv.  In the context of vector integration, it doesn't make sense to provide these.

### Test with sample data from Polyidus (hpv integration)

First, I tested ViFi with the sample data from Polyidus `https://www.pmgenomics.ca/hoffmanlab/proj/polyidus/polyidus-data-v1.tar.gz`.

I used these reads together with the referencences and HMMs provided with ViFi - 
`https://drive.google.com/open?id=0ByYcg0axX7udUDRxcTdZZkg0X1k`
`https://drive.google.com/open?id=0Bzp6XgpBhhghSTNMd3RWS2VsVXM`

```
bash 0_sample-data.sh
```
When running this, there was an error after bwa alignment:
```
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -t 1 -M /home/repo/data//hpv/hg19_hpv.fas /home/fastq/SiHa_R1.fastq.gz /home/fastq/SiHa_R2.fastq.gz
[main] Real time: 7.387 sec; CPU: 6.695 sec
Traceback (most recent call last):
  File "/home/scripts/get_trans_new.py", line 238, in <module>
    miscFile.write(b)
AttributeError: 'NoneType' object has no attribute 'write'
Prepared sequences for searching against HMMs: 0.010445s
Running HMMs
```

This looks to be just a typo in line 238:

```
if (miscFile is not None and len([read for read in q2aligns if read.is_unmapped]) + len([read for read in q1aligns if read.is_unmapped])) == 0:
   ^
```

I think should be:
```
if miscFile is not None and (len([read for read in q2aligns if read.is_unmapped]) + len([read for read in q1aligns if read.is_unmapped])) == 0:
														^
```
I cloned the ViFi repo and made this change.  In order to get the new version of ViFi to run inside the container, I needed to do:
```
export SINGULARITYENV_VIFI_DIR="/scratch1/sco305/intvi_simulation-experiments/ViFi"
```

This variable is passed into the singularity container, and ViFi looks for the variable `$VIFI_DIR` and runs scripts from `$VIFI_DIR/scripts`.


### Test with simulated data (test-easier)

I next tried to run ViFi with some simulated data.  This data was simulated (and analysed) by running `src/experiment0_prelim/1_virus.sh`.  The data was simulated from integrations of rep68 into chr1, so I concatenated rep68 together with hg19 (hg38 would be better but I leave that for later) and built bwa indices, and then ran ViFi (with HMMs disabled):

```
bash 1_sim-data.sh
```

However, there are errors after bwa alignment:
```
[main] Version: 0.7.17-r1188
[main] CMD: bwa mem -t 1 -M /home/repo/data//rep68/hg19_rep68.fas /usr/share/cond0.rep01.fq /usr/share/cond0.rep02.fq
[main] Real time: 1.879 sec; CPU: 1.600 sec
460 117 205
Traceback (most recent call last):
  File "/home/scripts/merge_viral_reads.py", line 128, in <module>
    scores = read_scores_file(args.reducedName[0])
  File "/home/scripts/merge_viral_reads.py", line 21, in read_scores_file
    input = open(hmm_file, 'r')  
IOError: [Errno 2] No such file or directory: 'tmp/temp/reduced.csv'
0
[Running BWA]: 1.530433
[Finished BWA]: 3.431194
```

### debugging from the ViFi repo

It's difficult to debug the code inside the container, because the Dockerfile for the ViFi container is not provided on Docker Hub, and singularity containers are read-only.  So, instead, clone the ViFi repo and run this code inside the container.  ViFi uses the environment variable `$VIFI_DIR` to find the scripts.

Use the script `2_debug-vifi-singularity.sh` to debug the ViFi scripts - this drops the user into a bash shell inside the singularity container, with all the bind-mounts for data created, and the `$VIFI_DIR` environment variable set

```
2_debug-vifi-singularity.sh

# once srun job has started
VIFI_DIR="/scratch1/sco305/intvi_simulation-experiments/ViFi"
CPUS="1"
READ1="cond0.rep01.fq"
READ2="cond0.rep02.fq"
python $VIFI_DIR/scripts/run_vifi.py -c ${CPUS} -f /usr/share/${READ1} -r /usr/share/${READ2} -v rep68 -o /home/output/ -d True
```

This reproduces the error from above, that looks to me to be related to the fact that we're not using HMMs, but the main ViFi script `run_vifi.py` expects some files to be present that are only created when using HMM's.  In particular, it looks like there are some reads identified based on alignments, some reads identified based on HMMs, and then these reads are merged.  When the merging occurs, ViFi expects the reads found from HMMs to be present in a file `tmp/temp/reduced.csv`, even if we've disabled HMMs with the command line argument `-d True`:

```
[Cluster and identify integration points]: 3.399760
Traceback (most recent call last):
  File "/scratch1/sco305/intvi_simulation-experiments/ViFi/scripts/merge_viral_reads.py", line 128, in <module>
    scores = read_scores_file(args.reducedName[0])
  File "/scratch1/sco305/intvi_simulation-experiments/ViFi/scripts/merge_viral_reads.py", line 21, in read_scores_file
    input = open(hmm_file, 'r')  
IOError: [Errno 2] No such file or directory: 'tmp/temp/reduced.csv'
0
[Finished cluster and identify integration points]: 5.144033
```

This is caused by ViFi running the script `merge_viral_reads.py` as a subprocess with the following arguments:
```
python /scratch1/sco305/intvi_simulation-experiments/ViFi/scripts/merge_viral_reads.py --unknown output.unknown.bam --trans output.trans.bam --reduced tmp/temp/reduced.csv --map tmp/temp/unmapped.map --output output.fixed.trans.bam
```
To produce the bam file `output.fixed.trans.bam`.  If we make this merging conditional on HMMs being used, we can just use the input file `output.trans.bam` to skip the next step, by replacing lines 155 and 156 in `run_vifi.py`:

```
os.system("python %s/scripts/merge_viral_reads.py --unknown %s.unknown.bam --trans %s.trans.bam --reduced tmp/temp/reduced.csv --map tmp/temp/unmapped.map --output %s.fixed.trans.bam" % (vifi_dir, options.prefix, options.prefix, options.prefix))
os.system("samtools sort -m 2G -@ %d %s.fixed.trans.bam > %s.fixed.trans.cs.bam" % (options.cpus, options.prefix, options.prefix))
```
with the conditional statement
```
if options.disable_hmms is False:
	os.system("python %s/scripts/merge_viral_reads.py --unknown %s.unknown.bam --trans %s.trans.bam --reduced tmp/temp/reduced.csv --map tmp/temp/unmapped.map --output %s.fixed.trans.bam" % (vifi_dir, options.prefix, options.prefix, options.prefix))
	os.system("samtools sort -m 2G -@ %d %s.fixed.trans.bam > %s.fixed.trans.cs.bam" % (options.cpus, options.prefix, options.prefix))
else:        
	os.system("samtools sort -m 2G -@ %d %s.trans.bam > %s.fixed.trans.cs.bam" % (options.cpus, options.prefix, options.prefix))
```

Now ViFi runs without any errors.



