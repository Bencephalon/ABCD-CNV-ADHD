import pandas as pd
import os

#### Task 1: Filter out CNVs with more than 30 CNVs > 100 kb ####
wd="/data/NCR_SBRB/jungbt/CNV/ABCD"
cnv_file=f"{wd}/results/penncnv/merged/ABCD.rawcnv"
#Load PennCNV raw output
cnv_df = pd.read_csv(cnv_file,delim_whitespace=True,header=None,
                    names=["coordcnv","numsnp","length","cn","file","startsnp","endsnp","conf"])
cnv_df["id"] = cnv_df["file"].apply(lambda x: os.path.basename(x))
cnv_df["length"] = cnv_df["length"].apply(lambda x: int(x.split("=")[1].replace(",","")))
#Identify all CNVs larger than 100 kb
cnv_df = cnv_df.loc[cnv_df.length > 100000,:]
# Count number of large CNVs per individual
count_df = cnv_df["id"].value_counts()
# Identify individuals with > 30 large CNVs
count_df = count_df.loc[count_df > 30]
# Save individuals with high CNV counts
count_df.to_csv(f"{wd}/results/penncnv/merged/high_cnv_counts.csv")

#Set QC thresholds
qc_dict = {"LRR_SD":0.35,
           "BAF_drift":0.01,
           "abs_WF":0.05,}

#Load PennCNV QC values
qc_file = f"{wd}/results/penncnv/merged/ABCD.filtered.qcsum"
df = pd.read_csv(qc_file,delim_whitespace=True)
df["File"] = df["File"].apply(lambda x: os.path.basename(x))
df["abs_WF"] = df["WF"].apply(lambda x: abs(x))

pass_df = df.copy()

#Remove high count subjects
print(f"{count_df.shape[0]} subjects have > 30 CNVs larger than 100 kb")
pass_df = pass_df.loc[~pass_df.File.isin(count_df.index),:]

# Identify individuals that fail LRR SD QC
qc_crit = "LRR_SD"
qc_filter = df[qc_crit] > qc_dict[qc_crit]
print(f"{sum(qc_filter)} subjects failed {qc_crit}")
pass_df = pass_df.loc[pass_df[qc_crit] <= qc_dict[qc_crit],:]

# Identify individuals that fail BAF drift QC
qc_crit = "BAF_drift"
qc_filter = df[qc_crit] > qc_dict[qc_crit]
print(f"{sum(qc_filter)} subjects failed {qc_crit}")
pass_df = pass_df.loc[pass_df[qc_crit] <= qc_dict[qc_crit],:]

# Identify individuals that fail GC waviness factor QC
qc_crit = "abs_WF"
qc_filter = df[qc_crit] > qc_dict[qc_crit]
print(f"{sum(qc_filter)} subjects failed {qc_crit}")
pass_df = pass_df.loc[pass_df[qc_crit] <= qc_dict[qc_crit],:]

# Identify all individuals that fail QC
failed_df  = df.drop(pass_df.index,axis=0)
print(f"{failed_df.shape[0]} subjects failed combined QC")
print(failed_df)

#Save failed subjects
failed_df.File.apply(lambda x: "_".join(x.split("_")[-2:]).replace(".txt","")).to_csv(f"{wd}/results/penncnv/merged/fail_qc.csv",header=False,index=False)

#Save passed subjects
pass_df.File.apply(lambda x: "_".join(x.split("_")[-2:]).replace(".txt","")).to_csv(f"{wd}/results/penncnv/merged/pass_qc.csv",header=False,index=False)

def reformat_penncnv(PN_df,prefix='ABCD_'):
    """
    Reformat PennCNV style inputs
    """
    PN_df.columns = ["coordcnv","numsnp","length","cn","file","startsnp","endsnp","conf"]
    PN_df["chr"]      = PN_df["coordcnv"].apply(lambda x: x.split(":")[0])
    PN_df["start"] = PN_df["coordcnv"].apply(lambda x: x.split(":")[1].split("-")[0]).astype(int)
    PN_df["end"]  = PN_df["coordcnv"].apply(lambda x: x.split(":")[1].split("-")[1]).astype(int)
    PN_df["numsnp"]     = PN_df["numsnp"].apply(lambda x: x.split("=")[1]).astype(int)
    PN_df["length"]     = PN_df["length"].apply(lambda x: x.split("=")[1].replace(",","")).astype(int)
    PN_df = PN_df.loc[PN_df["length"] > 1,:]
    PN_df["cn"]         = PN_df["cn"].apply(lambda x: x.split("=")[1]).astype(int)
    PN_df["ID"]         = PN_df["file"].apply(lambda x: os.path.basename(os.path.splitext(x)[0]).replace(prefix,""))
    PN_df["startsnp"]   = PN_df["startsnp"].apply(lambda x: x.split("=")[1])
    PN_df["endsnp"]     = PN_df["endsnp"].apply(lambda x: x.split("=")[1])
    PN_df["conf"]       = PN_df["conf"].apply(lambda x: x.split("=")[1]).astype(float)
    return PN_df


#Convert to a CSV
cnv_df = pd.read_csv(cnv_file,delim_whitespace=True,header=None,
                    names=["coordcnv","numsnp","length","cn","file","startsnp","endsnp","conf"])
cnv_df = reformat_penncnv(cnv_df)
cnv_df.ID = cnv_df.ID.apply(lambda x: "_".join(x.split("_")[-2:]))
# Remove subjects that failed QC
cnv_df = cnv_df.loc[~cnv_df.file.apply(lambda x: os.path.basename(x)).isin(failed_df.File),:]
#Save CNV file
cnv_df.to_csv(f"{wd}/results/penncnv/merged/penncnv_CNVs.csv")
