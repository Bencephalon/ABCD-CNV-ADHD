import pandas as pd
import numpy as np
import os
from tqdm import tqdm
from glob import glob
import sys

#Define paths
wd = "/data/NCR_SBRB/jungbt/CNV/ABCD"
hg19_region_dir = "/data/NCR_SBRB/jungbt/CNV/hg19_regions"

#Determine which pipeline is being run
pipeline = "consensus"
methods = ["PN","QT"]

# Make output directories
outdir = f"{wd}/results/{pipeline}"
os.makedirs(outdir,exist_ok=True)

# Set random seed
seed = 2514
np.random.seed(seed)


def calc_overlap(start1,stop1,start2,stop2):
    """
    Calculate the base pair overlap of two sequences on the same chr
    """
    #Calculate the length of seq 1
    var_length = stop1 - start1
    #Calculate the overlap
    var_overlap = range(max(start1,start2),min(stop1,stop2))
    return len(var_overlap)/var_length

def filter(df1,df2,thresh = 0.5):
    """
    Remove entries from df1 that don't overlap with entries in df2
    to a proportion of thresh
    """
    #Loop through each chromosome
    for chr in tqdm(df1["chr"].unique()):
        #Limit CNVs to those on a given chromosome
        chr_df1 = df1.loc[df1["chr"] == chr,:].sort_values(by="start")
        chr_df2 = df2.loc[df2["chr"] == chr,:].sort_values(by="start")
        #Loop through each CNV in df1
        for i, cnv in chr_df1.iterrows():
            #Loop through each CNV in df2
            for j, region in chr_df2.iterrows():
                #Check for overlap between each set of CNVs
                olap = calc_overlap(cnv["start"],cnv["end"],region["start"],region["end"])
                # If the olap passes a certain threshold, remove the CNV
                if olap > thresh:
                    # print(f"Remove {cnv['coordcnv']}")
                    df1 = df1.drop(i,axis = 0)
                    break
                elif region.start > cnv.end:
                    break
    return df1

def filter_individual(df1,df2,thresh = 0.7,behavior="keep"):
    """
    Remove entries from df1 that don't overlap with entries in df2
    to a proportion of thresh. Used for checking consensus.
    """
    #Specify output df
    out_df = pd.DataFrame(columns=df1.columns)
    #Loop through each individual
    for id in tqdm(df1["ID"].unique()):
        #Loop through each chromosome
        for chr in df1.loc[df1["ID"] == id,"chr"].unique():
            #Limit CNVs to a given chromosome
            chr_df1 = df1.loc[np.logical_and(df1["chr"] == chr,df1["ID"] == id),:].sort_values(by="start")
            chr_df2 = df2.loc[np.logical_and(df2["chr"] == chr,df2["ID"] == id),:].sort_values(by="start")
            #Loop through CNVs in df1
            for i, cnv in chr_df1.iterrows():
                #Loop through CNVs in df2
                for i, region in chr_df2.loc[chr_df2["cn"] == cnv["cn"],:].iterrows():
                    #Calculate overlap between each set of CNVs
                    olap = calc_overlap(cnv["start"],cnv["end"],region["start"],region["end"])
                    #If there is sufficient overlap
                    if olap >= thresh:
                        #If specified to keep olapping CNVs, keep it
                        if behavior == "keep":
                            out_df.loc[i,:] = cnv
                            break
                        #If specified to discard olapping CNVs, discard it
                        elif behavior == "discard":
                            break
                    elif region.start > cnv.end:
                        if behavior == "discard":
                            out_df.loc[i,:] = cnv
                else:
                    # If no CNVs match a given CNV, and instructed to discard
                    # overlapping CNVs, then keep the CNV
                    if behavior == "discard":
                        out_df.loc[i,:] = cnv
    return out_df

def summary(df,step):
    """
    Print a summary of the remaining CNVs
    """
    print(step)
    print(f"{df.shape[0]} CNVs")
    print(f"{sum(df['cn'] > 2)} duplications")
    print(f"{sum(df['cn'] < 2)} deletions")
    print(f"{len(df['ID'].unique())} subjects")
    print("-------------------")

def merge_cnvs(df,thresh = 0.2):
    """
    Combine two CNVs into one if the number of basepairs between them divided by
    the total length of the two CNVs is less than [thresh]
    """
    #Create output df
    new_df = pd.DataFrame(columns = df.columns)
    i = 0
    # for sub in tqdm(df.ID.unique(),total=len(df.ID.unique())):
    for sub in df.ID.unique():
        #Limit CNVs to a given subject
        tmp_df = df.loc[df.ID == sub,:]
        for chr in tmp_df.chr.unique():
            #Limit CNVs to a given chromosome, sorted by starting bp
            tmp_df2 = tmp_df.loc[tmp_df.chr == chr,:].sort_values(by="start")
            #If there's only 1 CNV on a given chromosome in a subject, there' no
            #chance of merging. Just add it
            if tmp_df2.shape[0] == 1:
                new_df.loc[i,:] = tmp_df2.iloc[0,:]
                new_df.loc[i,"conf"] = 0
                i += 1
            else:
                # Loop through each CNV
                for j, cnv in tmp_df2.iterrows():
                    # If the CNV isn't the first CNV, check for merging with the previous one
                    if not j == tmp_df2.index[0]:
                        # If the CNV copy number matches, then a match may be possible
                        if cnv["cn"] == new_df.loc[i-1,"cn"]:
                            #Calculate the ratio between the total length and the gap
                            gap = (cnv["start"] - new_df.loc[i-1,"end"])/(cnv["end"] - new_df.loc[i-1,"start"])
                            # If the gap ratio is smaller than the threshold
                            if gap <= thresh:
                                #Merge the CNVs
                                print(f"Merging {new_df.loc[i-1,'coordcnv']} with {cnv['coordcnv']} in {sub}")
                                #Replace the endsnp and end bp with the new CNV boundary
                                new_df.loc[i-1,"end"]      = cnv["end"]
                                new_df.loc[i-1,"endsnp"]   = cnv["endsnp"]
                                #Recode the coordcnv
                                new_df.loc[i-1,"coordcnv"] = f"{new_df.loc[i-1,'chr']}:{new_df.loc[i-1,'start']}-{new_df.loc[i-1,'end']}"
                                #Add together the snp counts
                                new_df.loc[i-1,"numsnp"]   = new_df.loc[i-1,"numsnp"] + cnv["numsnp"]
                                #Calculate the new length
                                new_df.loc[i-1,"length"]   = cnv["end"] - new_df.loc[i-1,"start"]
                                #Indicate that the CNV is merged
                                new_df.loc[i-1,"conf"]     -= 1
                                continue
                    #If none of the merging criterea are met, add the unmerged CNV to the final list
                    new_df.loc[i,:] = tmp_df2.iloc[0,:]
                    new_df.loc[i,"conf"] = 0
                    i += 1
    #Return the newly merged
    return new_df

#Find subjects that have a valid PennCNV input file
pn_files = glob(f"{wd}/data/penncnv/*.txt")
pn_files = ["_".join(os.path.basename(x).replace(".txt","").split("_")[-2:]) for x in pn_files]

#Load subject demographic and diagnosis sheet
sub_df = pd.read_csv(f"{wd}/data/assessments/subject_data.csv",index_col=0)

#Make a single variable for BATCH and METHOD (will save time later)
sub_df["BATCH_METHOD"] = sub_df.BATCH + "_" + sub_df.METHOD

# Remove subjects that are lacking either a demographic entry or PennCNV file
sub_df = sub_df.loc[sub_df.index.isin(pn_files),:]

#Load data from QuantiSNP and PennCNV with valid demographic data
if "PN" in methods:
    pn_df = pd.read_csv(f"{wd}/results/penncnv/merged/penncnv_CNVs.csv",index_col=0)
    pn_df  = pn_df.loc[pn_df.ID.isin(sub_df.index),:]

if "QT" in methods:
    qt_df = pd.read_csv(f"{wd}/results/quantisnp/cleaned/quantisnp_CNVs.csv",index_col=0)
    qt_df  = qt_df.loc[qt_df.ID.isin(sub_df.index),:]

#Load failed subjects for both pipelines to make sure subject selection is consistent
fail = []
fail += list(np.loadtxt(f"{wd}/results/penncnv/merged/fail_qc.csv",dtype=str))
fail += list(np.loadtxt(f"{wd}/results/quantisnp/cleaned/fail_qc.csv",dtype=str))

#Remove failed subjects
fail = set(fail)
sub_df = sub_df.loc[~sub_df.index.isin(fail),:]

if "PN" in methods:
    #Remove failed subjects
    pn_df  = pn_df.loc[~pn_df.ID.isin(fail),:]
    #Remove X chromosome
    pn_df  = pn_df.loc[~pn_df.chr.isin(["chrX"]),:]
    #Remove CNVs smaller than 20 kb (speeds things up)
    pn_df = pn_df.loc[pn_df["length"] >= 20000,:]
    # Print summary of CNVs
    summary(pn_df,"PennCNV CNVs:")

if "QT" in methods:
    #Remove failed subjects
    qt_df  = qt_df.loc[~qt_df.ID.isin(fail),:]
    #Remove X chromosome
    qt_df  = qt_df.loc[~qt_df.chr.isin(["chrX"]),:]
    #Remove CNVs smaller than 5 kb (speeds things up)
    qt_df = qt_df.loc[qt_df["length"] >= 20000,:]
    # Print summary of CNVs
    summary(qt_df,"QuantiSNP CNVs:")

z = 1
print(f"Step {z}: Remove WB Samples due to evidence of CNV differences")
z += 1
#NOTE: REMOVE ALL SUBJECTS FROM WHOLE-BLOOD SAMPLE
sub_df = sub_df.loc[~(sub_df.METHOD == "WB"),:]

#Incorporate QC Metrics
#Read in QC values
file_qc = f"{wd}/results/penncnv/merged/ABCD.filtered.qcsum"
qc_df = pd.read_csv(file_qc,sep="\t")
qc_df.index = qc_df.File.apply(lambda x: "_".join(os.path.basename(x).replace(".txt","").split("_")[-2:]))
#Assign QC values to subject df
sub_df["LRR_SD"]    = qc_df.LRR_SD
sub_df["BAF_SD"]    = qc_df.BAF_SD
sub_df["BAF_drift"] = qc_df.BAF_drift
sub_df["GCWF_abs"]  = abs(qc_df.WF)

#Save Subject DF
sub_df.to_csv(f"{outdir}/subject_info.csv")

#Threshold CNVs by confidence score
print(f"Step {z}: Threshld CNVs Based on Confidence")
z += 1
t = 20
qt_df = qt_df.loc[qt_df["conf"] >= t,:]
summary(qt_df,f"QuantiSNP CNVs (Conf >= {t}):")
pn_df = pn_df.loc[pn_df["conf"] >= t,:]
summary(pn_df,f"PennCNV CNVs (Conf >= {t}):")

print(f"Step {z}: Merge Adjacent CNVs")
z += 1

#Merge adjacent CNVs
if "PN" in methods:
    pn_df = merge_cnvs(pn_df)
    summary(pn_df,"PennCNV CNVs (merged):")

if "QT" in methods:
    qt_df = merge_cnvs(qt_df)
    summary(qt_df,"QuantiSNP CNVs (merged):")

#Check for consensus by removing CNVs that are not present in both quantiSNP and PennCNV
print(f"Step {z}: Check for Consensus")
df = filter_individual(pn_df,qt_df)
summary(df,f"Step {z}: Check for Consensus")
z += 1

print(f"Step {z}: Filter CNVs by Length and NumSNPS")
#Filter by length and numsnp
snp_thresh = 10
length_thresh = 50000
#Filter CNVs
df = df.loc[df["numsnp"] >= snp_thresh,:]
df = df.loc[df["length"] >= length_thresh,:]
summary(df,f"Step {z}: Filter CNVs by Length and NumSNPS")
z += 1

print(f"Step {z}: Remove CNV calls from telomeres and centromere")
telocentromeric_regions = pd.read_csv(f"{hg19_region_dir}/telomere_and_centromere_regions.csv")
#If a region is less than half of the length as our length threshold, we don't have to check it
telocentromeric_regions["length"] = telocentromeric_regions["chromEnd"] - telocentromeric_regions["chromStart"]
telocentromeric_regions.loc[telocentromeric_regions["length"] > length_thresh/2,:]
telocentromeric_regions.columns = ["chr","start","end","length"]
#Filter our CNVs by telomeric and centromeric regions
df = filter(df,telocentromeric_regions)
summary(df,f"Step {z}: Remove CNV calls from certain regions")
z += 1

print(f"Step {z}: Remove CNV calls from segmental duplication regions and MHC")
segmental_duplication_regions = pd.read_csv(f"{hg19_region_dir}/segmental_duplication_regions.csv")
#If a region is less than half of the length as our length threshold, we don't have to check it
segmental_duplication_regions["length"] = segmental_duplication_regions["chromEnd"] - segmental_duplication_regions["chromStart"]
# print(segmental_duplication_regions.shape[0])
segmental_duplication_regions.loc[segmental_duplication_regions["length"] > length_thresh/2,:]
# print(segmental_duplication_regions.shape[0])

MHC_region = pd.read_csv(f"{hg19_region_dir}/segmental_duplication_regions.csv")
exclude_50_df = pd.concat([segmental_duplication_regions,MHC_region])
exclude_50_df.columns = ["chr","start","end","length"]
df = filter(df,exclude_50_df)
summary(df,f"Step {z}: Remove CNV calls from segmental duplication regions and MHC")
z += 1

#Check the DGV for overlapping variants
print(f"Step {z}: Remove CNV calls overlapping common variants")
common = pd.read_csv(f"{hg19_region_dir}/common_variants_hg19.csv",index_col = 0)
#If a region is less than half of the length as our length threshold, we don't have to check it
common["length"] = common["chromEnd"] - common["chromStart"]
common.loc[common["length"] > length_thresh/2,:]
print(common)
common.columns = ["chr","start","end","cnv_type","size_cat","length"]
#Filter out duplicates
common_dup = common.loc[common.cnv_type.isin(["dup","both"])]
df_dup     = df.loc[df.cn > 2,]
df_dup = filter(df_dup,common_dup)
#Filter out deletions
common_del = common.loc[common.cnv_type.isin(["del","both"])]
df_del     = df.loc[df.cn < 2,]
df_del = filter(df_del,common_del)

df = pd.concat([df_del,df_dup],axis=0)

summary(df,f"Step {z}: Remove CNV calls overlapping common variants")
z += 1
####################### 4/12/2022 #######################

################################# SAVING DATA #################################

#Save CNV Calls in a CSV
print(f"Step {z}: Save CNV calls")
df.to_csv(f"{outdir}/processed_cnvs.csv",index=False)
