import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
import os
import seaborn as sns

# Paths
wd="/data/NCR_SBRB/jungbt/CNV/ABCD"

#Specify pipeline
pipeline = "consensus"
input = pipeline

#Specify number of ancestral populations (based on CV from ADMXITURE)
K=6

#Load subject information
sub_df = pd.read_csv(f"{wd}/results/{input}/subject_info.csv")
sub_df.index = sub_df.ID
print(sub_df)
print(sub_df.index)

#Load PCs
pc_df = pd.read_csv(f"{wd}/results/PCA/whole_sample/pca_analysis.eigenvec",sep = "\t",index_col=1)
pc_df_wnh = pd.read_csv(f"{wd}/results/PCA/WNH/pca_analysis.WNH.eigenvec",sep = "\t",index_col=1)

#remove subjects where pcs were not generated
sub_df = sub_df.loc[sub_df.index.isin(pc_df.index),:]

pc_df = pc_df.loc[sub_df.index,]
pc_df = pc_df.loc[~pc_df.index.duplicated(),:]

pc_df_wnh = pc_df_wnh.loc[~pc_df_wnh.index.duplicated(),:]


# Add PCs to subject df
sub_df[[f"PC{i}" for i in range(1,11)]] = pc_df[[f"PC{i}" for i in range(1,11)]]
# Add PCs for WNH individuals
pc_df_wnh.columns = [f"{x}.WNH" for x in pc_df_wnh.columns]
sub_df.loc[pc_df_wnh.index,[f"PC{i}.WNH" for i in range(1,11)]] = pc_df_wnh.loc[:,[f"PC{i}.WNH" for i in range(1,11)]]


#Load ADMIXTURE Ancestry
ancestry_df = pd.read_csv(f"{wd}/results/PCA/whole_sample/snp_for_ancestry.{K}.Q",header=None,sep=" ")
ancestry_df.index = pd.read_csv(f"{wd}/results/PCA/whole_sample/snp_for_ancestry.fam",header=None,sep="\t")[1]
ancestry_df = ancestry_df.loc[sub_df.index,]
ancestry_df = ancestry_df.loc[~ancestry_df.index.duplicated(),:]

#Define ancestry based on probs
ancestry_df["Ancestry"] = ancestry_df.idxmax(axis=1)

#Add ancestry to subject df
sub_df["Ancestry"] = ancestry_df["Ancestry"]


#Remove NAs
sub_df = sub_df.loc[~sub_df[[f"PC{i}" for i in range(1,11)] + ["Ancestry","sex","BATCH_METHOD"]].isnull().any(axis=1),:]

sub_df["Ancestry"] = "Ancestry_" + sub_df["Ancestry"].astype(int).astype(str)
print(sub_df.Ancestry.value_counts())

#Save all subject info
sub_df.to_csv(f"{wd}/results/{pipeline}/phenotypes_all.csv",index=False)

#Save revised subject info
sub_df.to_csv(f"{wd}/results/{pipeline}/phenotypes.csv",index=False)

#Generate plot of PCs vs Ancestry
for metric in ["Ancestry"]:
    plt.figure(figsize=[8,8])
    n=3
    #Plot Ancestry by PCs
    #PC1 x PC2
    plt.subplot(n,n,1)
    sns.scatterplot(x="PC1",y="PC2",hue=metric,data=sub_df,legend=False,size=0.1)
    plt.xticks([], [])
    plt.yticks([], [])

    #PC1 x PC3
    plt.subplot(n,n,2)
    sns.scatterplot(x="PC1",y="PC3",hue=metric,data=sub_df,legend=False,size=0.1)
    plt.xticks([], [])
    plt.yticks([], [])

    #PC1 x PC4
    plt.subplot(n,n,3)
    sns.scatterplot(x="PC1",y="PC4",hue=metric,data=sub_df,legend=False,size=0.1)
    plt.xticks([], [])
    plt.yticks([], [])

    #PC2 x PC3
    plt.subplot(n,n,5)
    sns.scatterplot(x="PC2",y="PC3",hue=metric,data=sub_df,legend=False,size=0.1)
    plt.xticks([], [])
    plt.yticks([], [])

    #PC2 x PC4
    plt.subplot(n,n,6)
    sns.scatterplot(x="PC2",y="PC4",hue=metric,data=sub_df,legend=False,size=0.1)
    plt.xticks([], [])
    plt.yticks([], [])

    #PC3 x PC4
    plt.subplot(n,n,9)
    sns.scatterplot(x="PC3",y="PC4",hue=metric,data=sub_df,legend=False,size=0.1)
    plt.xticks([], [])
    plt.yticks([], [])

    #PC3 x PC4
    plt.subplot(n,n,7)
    sns.scatterplot(x="PC3",y="PC4",hue=metric,data=sub_df,legend="brief",size=0,alpha=0)
    plt.axis("off")

    plt.tight_layout()
    plt.savefig(f"{wd}/results/{pipeline}/{metric}.png")
    plt.clf()
