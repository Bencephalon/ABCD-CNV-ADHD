#!/bin/bash
#Use swarm to run quantisnp

wd=/data/NCR_SBRB/jungbt/CNV/ABCD
cd ${wd}

mkdir logs/ABCD_quantisnp/

 swarm -g 2 -t 1 \
      --time 12:00:00 \
      --merge-output --logdir logs/ABCD_quantisnp/ \
      --partition norm --file run_quantisnp_ABCD.swarm
