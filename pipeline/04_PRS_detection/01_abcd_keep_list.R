library(tidyverse)
library(glue)
library(stringr)
source("scripts/99_util/util_load_cnvs.R")

#Load PLINK FAM file
plink_df <- read.table("data/PLINK/ABCD_release_3.0_QCed.filtered.fam",sep=" ")
colnames(plink_df) <- c("FID","IID","PID","MID","SEX","PHENO")
rownames(plink_df) <- plink_df$IID
print(dim(plink_df))

load_abcd_samples_init <- function(pipeline,pipeline_filter,opf = TRUE) {
  #Load subject information
  sub_df <- read_csv(file.path("results",pipeline,"subject_info.csv"),show_col_types = FALSE)
  sub_df <- sub_df[!startsWith(sub_df$ever_ADHD,"EXCLUDE_461"),]

  #NEW 07/07/2022: If a base file already exists, then take those individuals.
  # Useful for rerunning pipeline with same individuals
  if (file.exists(glue("results/consensus/data_files/{pipeline_filter}.rds"))) {
    print("Preselecting subjects...")
    keep_df <- readRDS(glue("results/consensus/data_files/{pipeline_filter}.rds"))
    #Remove all subjects that have a family member in keep_df, but aren't themselves in keep_df... Buh that's messy code
    sub_df <- sub_df[!(sub_df$rel_family_id %in% unique(keep_df$fam)) | sub_df$ID %in% rownames(keep_df),]
  }
  #Define various variables
  sub_df$ever_ADHD     <- factor(sub_df$ever_ADHD,levels=c("non-ADHD","ADHD"))
  sub_df$sex           <- factor(sub_df$sex)
  sub_df$BATCH_METHOD  <- factor(sub_df$BATCH_METHOD)
  sub_df$bipolar       <- factor(sub_df$ever_bipolar   == 1,levels=c(FALSE,TRUE))
  sub_df$psychosis     <- factor(sub_df$ever_psychosis == 1,levels=c(FALSE,TRUE))
  sub_df$bpd_psych     <- factor(sub_df$psychosis == TRUE | sub_df$bipolar == TRUE,levels=c(FALSE,TRUE))

  #Add variable for family
  sub_df$fam <- sub_df$rel_family_id
  #Remove subs with no stated family
  sub_df     <- sub_df[!is.na(sub_df$fam),]


  sub_df <- data.frame(sub_df)
  rownames(sub_df) <- sub_df$ID

  if(opf == TRUE) {
    #Select one member per family
    sub_df <- one_per_family_sample_targeted(sub_df)
  }

  for (status in c("ADHD","non-ADHD")) {
    print(status)
    print(glue("\tN = {sum(sub_df$ever_ADHD == status,na.rm = TRUE)}"))
    print(glue("\t%_Female = {sum(sub_df$ever_ADHD == status & sub_df$sex == 'F',na.rm = TRUE)/sum(sub_df$ever_ADHD == status,na.rm = TRUE)*100}"))
  }

  return(sub_df)
}

#Load in a list of subjects
keep_df <- load_abcd_samples_init("consensus","bpd_psych")

#Select only PLINK rows that correspond to our subset
plink_df <- plink_df[plink_df$IID %in% rownames(keep_df),]

#Save our analysis subset for PCA
write.table(plink_df[c("FID","IID")],"data/PLINK/keep_subjects_pca.txt",sep="\t",col.names=FALSE,row.names=FALSE,quote = FALSE)
