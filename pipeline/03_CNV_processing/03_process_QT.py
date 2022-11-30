#Goal: Convert QuantiSNP format to PennCNV format
#Imports
import os
import pandas as pd
import numpy as np
from glob import glob
from tqdm import tqdm
import shutil

#Specify paths
wd              = "/data/NCR_SBRB/jungbt/CNV/ABCD"

#Specify output directory, delete it if it already exists, and recreate it
QT_output_path  = f"{wd}/results/quantisnp"
shutil.rmtree(f"{QT_output_path}/cleaned/",ignore_errors=True)
os.makedirs(f"{QT_output_path}/cleaned/",exist_ok = True)

#Find QuantiSNP files
paths = glob(f"{QT_output_path}/*")
skip = 0
samples = []
# Combine all CNV calls into a single CNV file
for i, path in tqdm(enumerate(paths),total=len(paths)):
    sample = '_'.join(path.split("_")[-2:])
    samples.append(sample)
    #QC Header is incorrectly formatted. Repair it.
    if not os.path.exists(f"{path}/{sample}.qcnew"):
        try:
            with open(f"{path}/{sample}.qc","r+") as qcfile:
                qc = qcfile.readlines()
        except FileNotFoundError:
            print(f"{sample} output from QuantiSNP not found. Skipping.")
            skip += 1
            continue
        #Add gender to the header
        qc[0] = 'Sample ID\tChromosome\tOutlier Rate\tStd. Dev. LRR\tStd. Dev. BAF\tGender\n'
        #Remove gender from the first line
        qc[1] = qc[1].replace("Gender\t","")
        #Rewrite qc
        with open(f"{path}/{sample}.qcnew","w") as qcfile:
            for line in qc:
                qcfile.write(line)
    if i == 0:
        #Define merged Dataframes
        CNV_df  = pd.read_csv(f"{path}/{sample}.cnv",sep="\t")
        QC_df   = pd.read_csv(f"{path}/{sample}.qcnew",sep="\t")
    else:
        #Add CNVs to CNV file
        try:
            CNV_df  = CNV_df.append(pd.read_csv(f"{path}/{sample}.cnv",sep="\t"),ignore_index = True)
        except FileNotFoundError:
            print(f"{sample} output from QuantiSNP not found. Skipping.")
            skip += 1
            continue
        QC_df   = QC_df.append(pd.read_csv(f"{path}/{sample}.qcnew",sep="\t"),ignore_index = True)

print(f"{skip} subjects skipped")
print(f"{CNV_df.shape[0]} CNVs")
# Remove subjects with > 30 CNVs > 100 kb
cnv_len_df = CNV_df.copy()
# Identify CNVs larger than 100 kb
cnv_len_df = cnv_len_df.loc[cnv_len_df["Length (bp)"] > 100000,:]
# Count number of large CNVs per individual
count_df = cnv_len_df["Sample Name"].value_counts()
# Identify individuals with > 30 large CNVs
count_df = count_df.loc[count_df > 30]

#Save individuals with many large CNVs
count_df.to_csv(f"{wd}/results/quantisnp/cleaned/high_cnv_counts.csv")

#Define subjects that pass QC
pass_array = np.array(samples)
#Remove high count subjects
print(f"{count_df.shape[0]} subjects have > 30 CNVs larger than 100 kb")
pass_array = pass_array[~np.in1d(pass_array,count_df.index)]

#Define subjects that failed QC
failed_array = np.array(count_df.index)

#Save failed subjects
np.savetxt(f"{wd}/results/quantisnp/cleaned/fail_qc.csv",failed_array,fmt="%s")
#Save passed subjects
np.savetxt(f"{wd}/results/quantisnp/cleaned/pass_qc.csv",pass_array,fmt="%s")

#Convert CNVs to a consistent format for consensus matching with PennCNV
cnv_out_df = pd.DataFrame(columns = ['coordcnv', 'numsnp', 'length', 'cn', 'file', 'startsnp', 'endsnp', 'conf', 'chr', 'start', 'end', 'ID'])
cnv_out_df["chr"]       = [f"chr{i}" if i < 23 else "chrX" for i in CNV_df.Chromosome]
cnv_out_df["start"]     = CNV_df["Start Position (bp)"]
cnv_out_df["end"]       = CNV_df["End Position (bp)"]
cnv_out_df["coordcnv"]  = cnv_out_df.chr + ":" + cnv_out_df.start.astype(str) + "-" + cnv_out_df.end.astype(str)
cnv_out_df["numsnp"]    = CNV_df['No. Probes']
cnv_out_df["length"]    = CNV_df['Length (bp)']
cnv_out_df["cn"]        = CNV_df['Copy Number']
cnv_out_df["file"]      = CNV_df['Sample Name']
cnv_out_df["startsnp"]  = CNV_df['Start Probe ID']
cnv_out_df["endsnp"]    = CNV_df['End Probe ID']
cnv_out_df["conf"]      = CNV_df['Max. Log BF']
cnv_out_df["ID"]        = CNV_df['Sample Name']

#Exclude CNVs from subjects that failed QC
cnv_out_df = cnv_out_df.loc[~cnv_out_df.ID.isin(failed_array),:]

#Save CNVs
cnv_out_df.to_csv(f"{wd}/results/quantisnp/cleaned/quantisnp_CNVs.csv")
