#!/bin/bash

set -e

CONTAINER="rscripts"
VERSION="6"

if [ ! -e ${CONTAINER}_${VERSION}.sif ] ; then
	singularity pull docker://szsctt/${CONTAINER}:${VERSION}
fi

# run R notebook
singularity exec ${CONTAINER}_${VERSION}.sif Rscript -e 'rmarkdown::render("analysis/2021-11-29_simulation-figures.Rmd")'

