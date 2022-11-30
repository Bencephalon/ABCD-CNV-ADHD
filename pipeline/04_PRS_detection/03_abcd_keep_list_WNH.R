library(tidyverse)
library(glue)
library(stringr)
source("scripts/99_util/util_load_cnvs.R")

K = 6

#Load ADMIXTURE Ancestry
ancestry_df <- read.table(glue("results/PCA/whole_sample/snp_for_ancestry.{K}.Q"),header=FALSE,sep=" ")
fam_df <- read.table("results/PCA/whole_sample/snp_for_ancestry.fam",header=FALSE,sep="\t")
colnames(fam_df) <- c("FID","IID","PID","MID","SEX","PHENO")
rownames(ancestry_df) <- fam_df$IID
colnames(ancestry_df) <- paste("Ancestry",1:K,sep="_")

#Define ancestry based on probs
ancestry_df$Ancestry <- apply(ancestry_df,1,function(row) colnames(ancestry_df)[which.max(row)])

print(table(ancestry_df$Ancestry))

#Remove subjects from non-largest ancestry
largest_ancestry <- names(which.max( table(ancestry_df$Ancestry)))
fam_df <- fam_df[fam_df$IID %in% rownames(ancestry_df[ancestry_df$Ancestry == largest_ancestry,]),]

#Save white, non-Hispanic subjects
write.table(fam_df[c("FID","IID")],"data/PLINK/keep_subjects_pca.WNH.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote = FALSE)
