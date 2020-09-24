## Polyidus test scripts

Some scripts to test out [Polyidus](https://github.com/hoffmangroup/polyidus/tree/master/src).

`0_dependencies.sh`

Uses `mamba` to create a `conda` environment containing the dependencies required to run Polyidus.

`1_demo-data.sh`

Run Polyidus using the [demo data](https://www.pmgenomics.ca/hoffmanlab/proj/polyidus/polyidus-data-v1.tar.gz) provided by the developers, using the conda environment.

`2_run-test.sh`

Run Polyidus using some simulated data. Note that this requires the script `src/sim_tests/test-sim-analysis_combined.sh` to have already been run.

`3_local-debug.sh`

Run Polyidus on the same simulated data, using docker, on my local machine.  I made a docker container with the Polyidus requirements and scripts (docker://szsctt/polyidus).

`4_test-singularity.sh`

I made a docker container with Polyidus - test this out on the cluster using singularity.


