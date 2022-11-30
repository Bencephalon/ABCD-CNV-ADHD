import pandas as pd
import numpy as np
from tqdm import tqdm
from glob import glob
import os
import multiprocessing as mp
import time

start_time = time.time()

wd = "/data/NCR_SBRB/jungbt/CNV/ABCD"
outdir=f"{wd}/data/penncnv"

os.makedirs(outdir,exist_ok=True)

#Load list of SNPs to use for CNV detection
snplist = f"{wd}/data/genotype_QCed/plink.snplist"
snp_df = pd.read_csv(snplist,header=None,names=["SNP"])

def save_penncnv_file(id):
    """
    Reformat data and save it in a PennCNV-compatible format
    """
    #Check if File exists
    if not os.path.exists(f"{outdir}/{id}.txt"):
        #Create output df
        out_df = pd.DataFrame(index=range(lrr_df.shape[0]),columns=["Name","Chr","Position","B Allele Freq","Log R Ratio"])
        #Insert probe info
        out_df[["Name","Chr","Position"]] = probe_df[["dbSNP.RS.ID","Chromosome","Physical.Position"]]
        #Insert CNV metrics
        out_df["Log R Ratio"] = lrr_df[id]
        out_df["B Allele Freq"] = baf_df[id]
        #Save
        out_df.to_csv(f"{outdir}/{id}.txt",sep="\t",index=False)

#Loop through each batch of SNPs provided by the ABCD
for file in glob(f"{wd}/data/ABCD/BATCH_*.sample"):
    #Extract Batch Name
    batch = os.path.basename(file).replace(".sample","")
    print(batch)
    #Load Sample Info
    df       = pd.read_csv(file,sep="\t",header = None,
                            names=["samp","id"])
    #Generate unique ID
    df["unique_id"] = batch + "_" + df.samp + "_" + df.id
    #Load Probe Info
    probe_df = pd.read_csv(f"{wd}/data/ABCD/{batch}.probe.info",sep="\t")
    #Load LRR File
    lrr_df   = pd.read_csv(f"{wd}/data/ABCD/{batch}_lrr.txt",sep="\t",
                            header = None,names=df.unique_id)
    #Load BAF File
    baf_df   = pd.read_csv(f"{wd}/data/ABCD/{batch}_baf.txt",sep="\t",
                            header = None,names=df.unique_id)
    #Remove SNP data not in list of QC'd SNPs
    lrr_df   = lrr_df.loc[probe_df["dbSNP.RS.ID"].isin(snp_df.SNP),:].reset_index()
    baf_df   = baf_df.loc[probe_df["dbSNP.RS.ID"].isin(snp_df.SNP),:].reset_index()
    probe_df = probe_df.loc[probe_df["dbSNP.RS.ID"].isin(snp_df.SNP),:].reset_index()
    #Parallel Saving of PennCNV Files
    print("\tSaving PennCNV Files...")
    pool = mp.Pool(processes=16)
    results = pool.map(save_penncnv_file, df.unique_id)

print("--- %s seconds ---" % (time.time() - start_time))
