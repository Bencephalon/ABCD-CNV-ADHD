#Imports
import pandas as pd
from glob import glob
import os
from tqdm import tqdm
import multiprocessing as mp
import numpy as np
import shutil

#Input file paths
wd               = f'/data/NCR_SBRB/jungbt/CNV/ABCD'
qt_dir           = f'/data/NCR_SBRB/jungbt/CNV/packages/quantisnp'
qt_datadir       = f'{wd}/data/quantisnp'
sex_file         = f'{wd}/data/assessments/subject_data.csv'
QT_outdir        = f'{wd}/results/quantisnp'
out_swarmfile    = f'{wd}/run_quantisnp_ABCD.swarm'


os.makedirs(QT_outdir,exist_ok=True)
os.makedirs(qt_datadir,exist_ok=True)

#QuantiSNP params
MCRROOT		    = f"{qt_dir}/opt/MATLAB/MATLAB_Compiler_Runtime/v79" # path to MCR Run-Time Libraries
EMITERS         = 10				                                # number of EM iterations to use during training
LSETTING        = 2000000						                    # characteristic CNV length parameter
GCDIR           = f"{qt_dir}/gc_data_b37"                           # path to GC data files (contents of gc_data.zip)
SUBSAMPLELEVEL  = 1 							                    # number of sub-samples to use
PARAMSFILE      = f"{qt_dir}/config/params.dat"                     # path to parameters file
LEVELSFILE      = f"{qt_dir}/config/levels-hd.dat"                     # path to levels file
CHRRANGE        = "1:23"                                            # path to parameters file
CHRX            = 23                                                # which chromosome is X?
#Sample-specific variables: OUTDIR, ID, SEX, FILE

#Sex conversion dictionaries
sex_code   = {"M":"male","F":"female"}

# Open demographic file containing the sex and extract it
sex_df = pd.read_csv(sex_file,index_col=0)
sex_df = sex_df.loc[~sex_df.sex.isnull(),:]
sex_df.sex = sex_df.sex.apply(lambda x: sex_code[x])

#Open the samples that passed PennCNV QC
pass_qc = np.loadtxt(f"{wd}/results/penncnv/merged/pass_qc.csv",dtype=str)

#Get list of PennCNV data files
samples = glob(f"{wd}/data/penncnv/*.txt")
samples = [x for x in samples if "_".join(os.path.basename(x).split("_")[-2:]).replace(".txt","") in pass_qc]

#Generate swarm commands
swarm_lines = []
def process_cmds(sample):
    """
    Recode PennCNV data files for use in QuantiSNP.
    """
    # Identify the sample ID
    ID_OUT = os.path.splitext(os.path.basename(sample))[0]
    ID     = "_".join(ID_OUT.split("_")[-2:])
    # Define the input penncnv file
    FILE=f"{qt_datadir}/{ID_OUT}.txt"
    # Define the recoded file name output directory
    OUTDIR = f"{QT_outdir}/{ID_OUT}"
    # Check to see if the subject has a demographic entry
    # Skip if lacking
    try:
        SEX    = sex_df.loc[ID,"sex"]
    except KeyError:
        print(f"{ID} not found in data file. Skipping")
        return ""
    # Write out the QuantiSNP input file only if it doesn't exist
    if not os.path.exists(f"{wd}/data/quantisnp/{ID_OUT}.txt"):
        #Reformat PennCNV format
        df = pd.read_csv(sample,sep="\t",
                        dtype={'Name'    : str,
                               'Chr'          : str,
                               'Position'     : int,
                               'B Allele Freq': float,
                               'Log R Ratio'  : float})
        # Reformat file
        df = df[["Name","Chr","Position",f"Log R Ratio",f"B Allele Freq"]]
        df["Sample ID"] = ID
        #Save quantiSNP data file
        df.to_csv(f"{wd}/data/quantisnp/{ID_OUT}.txt",sep="\t",index=False)


    #Check if the QuantiSNP command has already been run by looking for an output file
    if os.path.exists(f"{OUTDIR}/{ID}.cnv"):
        return ""
    else:
        # Output the QuantiSNP command
        return f"mkdir {OUTDIR}; echo {ID}; {qt_dir}/linux64/quantisnp2 {MCRROOT} --chr {CHRRANGE} --outdir {OUTDIR} --sampleid {ID} --gender {SEX} --emiters {EMITERS} --lsetting {LSETTING} --gcdir {GCDIR} --config {PARAMSFILE} --levels {LEVELSFILE} --input-files {FILE} --chrX {CHRX} --doXcorrect --isaffy\n"

#Run process_cmds in parallel
pool = mp.Pool(64)
swarm_lines = pool.map(process_cmds, samples)
pool.close()
pool.join()
#Remove blanks
swarm_lines = [x for x in swarm_lines if not x == ""]

#Specify number of jobs to run
jobs = min(250,len(swarm_lines))

#Make output directory for batch jobs
os.makedirs(f"{qt_datadir}/batches",exist_ok=True)
#Delete previous batch jobs
for file in glob(f"{qt_datadir}/batches/*.sh"):
    try:
        os.remove(file)
    except:
        print("Error while deleting file : ", file)

#Split the data files into equally sized jobs
sub_files = np.array_split(swarm_lines,jobs)

# Write QuantiSNP commands to run CNV detection in batches
with open(out_swarmfile, "w") as writefile:
    for i, cmds in enumerate(sub_files):
        with open(f"{qt_datadir}/batches/batch_{i}.sh","w") as batchfile:
            #Write out the pre-reqs for QuantiSNP
            batchfile.write("module load matlab;\n")
            batchfile.write("LD_LIBRARY_PATH=.:/data/NCR_SBRB/jungbt/CNV/packages/quantisnp/opt/MATLAB/MATLAB_Compiler_Runtime/v79/runtime/glnxa64:/data/NCR_SBRB/jungbt/CNV/packages/quantisnp/opt/MATLAB/MATLAB_Compiler_Runtime/v79/bin/glnxa64:/data/NCR_SBRB/jungbt/CNV/packages/quantisnp/opt/MATLAB/MATLAB_Compiler_Runtime/v79/sys/os/glnxa64:/data/NCR_SBRB/jungbt/CNV/packages/quantisnp/opt/MATLAB/MATLAB_Compiler_Runtime/v79/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/data/NCR_SBRB/jungbt/CNV/packages/quantisnp/opt/MATLAB/MATLAB_Compiler_Runtime/v79/sys/java/jre/glnxa64/jre/lib/amd64/server:/data/NCR_SBRB/jungbt/CNV/packages/quantisnp/opt/MATLAB/MATLAB_Compiler_Runtime/v79/sys/java/jre/glnxa64/jre/lib/amd64/client:/data/NCR_SBRB/jungbt/CNV/packages/quantisnp/opt/MATLAB/MATLAB_Compiler_Runtime/v79/sys/java/jre/glnxa64/jre/lib/amd64;\n")
            batchfile.write("export LD_LIBRARY_PATH;\n")
            batchfile.write(f"XAPPLRESDIR={MCRROOT}/X11/app-defaults;\n")
            batchfile.write("export XAPPLRESDIR;\n")
            for cmd in cmds:
                #Write out the quantisnp commands
                batchfile.write(cmd)
        writefile.write(f"bash {qt_datadir}/batches/batch_{i}.sh\n")
