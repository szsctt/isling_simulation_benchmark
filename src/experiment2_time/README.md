## Scripts to test runtime

Simulate data with the simulation scripts, and then run various tools to test runtime and memory usage.  Run all jobs with exclsive use of one node, 20 cpus, 24hr max walltime, 128gb memory.  Use `/usr/bin/time` to get runtime and memory usage.  Use all tools inside singularity containers, pulled from my docker hub account:
 - isling: `docker://szsctt/isling:1`
 - ViFi: `docker://szsctt/vifi:1`
 - polyidus: `docker://szsctt/polyidus:3`
 - seeksv: `docker://szsctt/seeksv:1`
 
Top-level script is `run_time.sh`.  This uses the config file `../../config/experiment2_time/simulation.yml` to simulate data. 

It then runs `isling` and the snakemake workflow in `intvi_other-tools` to index the necessary references.  

It then runs a python script that does the timing - it first iterates over experiments in the config file, then samples simulated for that experiment.  For each sample, it runs each tool three times, and collects the results in a tab-separated format.


Note that many of the scripts assume that there is only one reference in the simulation config file!
