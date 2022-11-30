import pandas as pd
import numpy as np
from glob import glob


#Define number of swarm jobs to create
jobs = 100

#Define all file locations
wd="/data/NCR_SBRB/jungbt/CNV/ABCD"
hmm="/usr/local/apps/penncnv/1.0.4/lib/hhall.hmm"
gcmodel=f"{wd}/data/model/ABCD.hg19.gcmodel"
pfb=f"{wd}/data/model/ABCD.pfb"
swarm_file=f"{wd}/run_penncnv_ABCD.swarm"
out=f"{wd}/results/penncnv"
data_dir=f"{wd}/data/penncnv"

#Find PennCNV data files
files = np.array(glob(f"{data_dir}/*.txt"))
#Split PennCNV data files into equal-sized jobs
sub_files = np.array_split(files,jobs)

# Write PennCNV commands to run CNV detection in batches
with open(swarm_file, "w") as writefile:
    for i, sub_list in enumerate(sub_files):
        cmd = f"detect_cnv.pl -test -hmm {hmm} -gcmodel {gcmodel} -pfb {pfb} -confidence --minlength 20000 -log {out}/{i}_ABCD.log {' '.join(sub_list)} -out {out}/{i}_ABCD.rawcnv\n"
        writefile.write(cmd)
