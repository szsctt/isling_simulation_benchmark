#!/bin/bash

set -e

CONTAINER="rscripts"
VERSION="5"

if [ ! -e ${CONTAINER}_${VERSION}.sif ] ; then
	singularity pull docker://szsctt/${CONTAINER}:${VERSION}
fi

# run R notebook
singularity exec ${CONTAINER}_${VERSION}.sif Rscript -e 'rmarkdown::render("analysis/2021-06-15_new-isling-results.Rmd")'

