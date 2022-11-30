#!/bin/bash

wd=/data/NCR_SBRB/jungbt/CNV/ABCD

#Create output directory
mkdir -p ${wd}/results/penncnv/

#Run swarm batch commands
swarm -f ${wd}/run_penncnv_ABCD.swarm -m penncnv --time=4:00:00 --partition=norm,quick -g 5 -t 1
