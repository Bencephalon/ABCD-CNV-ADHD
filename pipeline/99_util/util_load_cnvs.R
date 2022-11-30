# Libaries
library(tidyverse)
library(glue)
library(stringr)

set.seed(2514)

setwd("/data/NCR_SBRB/jungbt/CNV/ABCD")

load_abcd_samples <- function(input,pipeline,pipeline_filter,sex = "all",dx="ksads") {
  #Load subject information
  sub_df <- read_csv(file.path("results",input,"phenotypes_all.csv"),show_col_types = FALSE)
  #Exclude plate 461
  sub_df <- sub_df[!startsWith(sub_df$ever_ADHD,"EXCLUDE_461"),]

  #07/07/2022: Take subjects from base file if it already exists. Don't pick randomly
  if (file.exists(glue("results/consensus/data_files/{pipeline_filter}.rds"))) {
    print("Preselecting subjects...")
    print(dim(sub_df))
    keep_df <- readRDS(glue("results/consensus/data_files/{pipeline_filter}.rds"))
    # If there's already a family member present, exclude all other members, otherwise keep all family members and we'll exclude others later
    sub_df <- sub_df[!(sub_df$rel_family_id %in% unique(keep_df$fam)) | sub_df$ID %in% rownames(keep_df),]
    # print(dim(sub_df))
  } else if (file.exists(glue("results/consensus/data_files/{pipeline_filter}_old.rds"))) {
    print("Preselecting subjects...")
    print(dim(sub_df))
    keep_df <- readRDS(glue("results/consensus/data_files/{pipeline_filter}_old.rds"))
    # If there's already a family member present, exclude all other members, otherwise keep all family members and we'll exclude others later
    sub_df <- sub_df[!(sub_df$rel_family_id %in% unique(keep_df$fam)) | sub_df$ID %in% rownames(keep_df),]
  } else if (file.exists(glue("data/PLINK/keep_subjects_pca.txt"))) {
    print("Preselecting subjects...")
    print(dim(sub_df))
    #Use the same subjects from the PCA analysis, which uses the same criteria
    keep_df <- read_table("data/PLINK/keep_subjects_pca.txt",col_names=FALSE,show_col_types = FALSE)
    sub_df <- sub_df[sub_df$ID %in% keep_df$X2,]
  }


  #Factorize several variables
  sub_df$ever_ADHD_str <- sub_df$ever_ADHD
  sub_df$ever_ADHD     <- factor(sub_df$ever_ADHD,levels=c("non-ADHD","ADHD"))
  sub_df$sex           <- factor(sub_df$sex)
  sub_df$BATCH_METHOD  <- factor(sub_df$BATCH_METHOD)
  sub_df$Ancestry      <- factor(sub_df$Ancestry)
  sub_df$externalizing <- factor(sub_df$externalizing,levels=c(FALSE,TRUE))
  sub_df$anxiety       <- factor(sub_df$anxiety,levels=c(FALSE,TRUE))
  sub_df$mood          <- factor(sub_df$mood,levels=c(FALSE,TRUE))
  sub_df$OCD           <- factor(sub_df$OCD,levels=c(FALSE,TRUE))
  sub_df$ODD           <- factor(sub_df$ODD == 1,levels=c(FALSE,TRUE))
  sub_df$CD            <- factor(sub_df$CD == 1,levels=c(FALSE,TRUE))
  sub_df$bipolar       <- factor(sub_df$ever_bipolar   == 1,levels=c(FALSE,TRUE))
  sub_df$psychosis     <- factor(sub_df$ever_psychosis == 1,levels=c(FALSE,TRUE))
  sub_df$bpd_psych     <- factor(sub_df$psychosis == TRUE | sub_df$bipolar == TRUE,levels=c(FALSE,TRUE))

  #Set CBCL thresholds
  cbcl_thresh_high = 65
  cbcl_thresh_low = 65

  #Define ADHD in CBCL
  sub_df$ADHD_CBCL <- sub_df$ADHD_CBCL_DSM5_Scale_t_max
  sub_df$ADHD_CBCL[sub_df$ADHD_CBCL >= cbcl_thresh_high] <- "ADHD"
  sub_df$ADHD_CBCL[sub_df$ADHD_CBCL < cbcl_thresh_low] <- "non-ADHD"
  sub_df$ADHD_CBCL    <- factor(sub_df$ADHD_CBCL,levels=c("non-ADHD","ADHD"))

  #Convert dx to KSAD-CBCL or CBCL if specified in pipeline
  if (dx == "cbcl") {
    print("Changing dx metric to CBCL")
    sub_df$ever_ADHD <- sub_df$ADHD_CBCL
    #Drop subjects without CBCl
    sub_df <- sub_df[!is.na(sub_df$ADHD_CBCL_DSM5_Scale_t_max),]
  } else if (dx == "inatt") {
    #Exclude hyperactive-only
    sub_df <- sub_df[!(sub_df$ever_ADHD == "ADHD" & sub_df$sx_inatt_present_max < 6),]
  } else if (dx == "HI") {
    sub_df <- sub_df[!(sub_df$ever_ADHD == "ADHD" & sub_df$sx_hi_present_max < 6),]
  } else if (dx == "ksad_cbcl") {
    print("Changing dx metric to agreement between KSADs and CBCL")
    print(head(sub_df$ADHD_CBCL))
    print(head(sub_df$ever_ADHD))
    sub_df <- sub_df[sub_df$ADHD_CBCL == sub_df$ever_ADHD,]
    sub_df <- sub_df[!is.na(sub_df$ADHD_CBCL_DSM5_Scale_t_max),]
    print("~~~~~~~~~~~~~~~~~")
    print(dim(sub_df))
    print("~~~~~~~~~~~~~~~~~")
  }

  #Add family
  sub_df$fam         <- sub_df$rel_family_id

  #Calculate the frequency of each ancestry to find the largest
  ancestry_freq <- sort(table(sub_df$Ancestry))

  # Limit sex if pipeline specifies
  if (sex == "female") {
    #Limit to females only
    sub_df <- sub_df[sub_df$sex == "F",]
  } else if (sex == "male") {
    #Limit to males only
    sub_df <- sub_df[sub_df$sex == "M",]
  }
  #Limit ancestry if pipeline specifies
  if (pipeline == "1st_eth") {
    #Limit to largest ancestry
    sub_df = sub_df[sub_df$Ancestry == names(ancestry_freq)[length(ancestry_freq)],]
    #Overwrite the PCs
    sub_df[paste("PC",1:10)] <- sub_df[glue("PC{1:10}.WNH")]
  }
  #Limit subjects based on other variables
  if (!(pipeline_filter %in% c("default","cbcl"))) {
    # Remove all subjects with 1 in the given criteria
    sub_df <- sub_df[sub_df[,pipeline_filter] == FALSE & !is.na(sub_df[,pipeline_filter]),]
  }
  #Format data frame
  sub_df <- data.frame(sub_df)
  rownames(sub_df) <- sub_df$ID
  #Add PRS Scores
  file.prs <- "/data/NCR_SBRB/jungbt/CNV/ABCD/PRS/results/ABCD_release_3.0_QCedfor_Ben.PRS_adhd2019.all_score"
  prs_measures <- c("Pt_5e.08","Pt_1e.05","Pt_0.0001","Pt_0.001","Pt_0.01","Pt_0.05","Pt_0.1","Pt_0.2","Pt_0.5","Pt_1")
  prs <- read.table(file.prs,sep="",header = TRUE,row.names = 2)
  #PRS analysis shows that P < 0.1 is most significantly associated with ADHD, so when using generic PRS, select that threshold
  sub_df$PRS <- prs[rownames(sub_df),"Pt_0.1"]
  #Add PRS measures
  sub_df[,prs_measures] <- prs[rownames(sub_df),prs_measures]

  #Set IQ as pea_wiscv_tss_max (IQ proxy)
  sub_df$IQ <- sub_df$pea_wiscv_tss_max
  sub_df$IQ_unscaled <- sub_df$pea_wiscv_tss_max

  #Select continuous variables to Z-score
  scale_cols <- c("IQ","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PRS",prs_measures)

  #Z-Scale our continuous data columns
  sub_df[scale_cols] <- scale(sub_df[scale_cols])

  return(sub_df)
}

load_abcd_burden <- function(input,pipeline,pipeline_filter,sex = "all",dx="ksads") {
  #Start by loading the relevant samples
  sub_df <- load_abcd_samples(input,pipeline,pipeline_filter,sex,dx=dx)
  #Load the pre-calculated CNV burden files
  burden_df <- read_csv(file.path("results",input,"analysis_tables","burden","burden_table.csv"),show_col_types = FALSE)
  colnames(burden_df)[1] <- "ID"
  #Recode each burden variable to include the type of CNV
  for (cnv_type in c("del","dup","both")) {
    burden_df <- burden_df[burden_df$ID %in% sub_df$ID,]
    burden_df[[glue("Binary_50kb_{cnv_type}")]]  <- factor(burden_df[[glue("Binary_50kb_{cnv_type}")]],levels=c(FALSE,TRUE))
    burden_df[[glue("Binary_100kb_{cnv_type}")]] <- factor(burden_df[[glue("Binary_100kb_{cnv_type}")]],levels=c(FALSE,TRUE))
    burden_df[[glue("Binary_200kb_{cnv_type}")]] <- factor(burden_df[[glue("Binary_200kb_{cnv_type}")]],levels=c(FALSE,TRUE))
  }
  #Limit output data frame to only a subset of columns
  prs_measures       <- c("Pt_5e.08","Pt_1e.05","Pt_0.0001","Pt_0.001","Pt_0.01","Pt_0.05","Pt_0.1","Pt_0.2","Pt_0.5","Pt_1")
  keep_cols.basic    <- c("ID","interview_age_max","BATCH_METHOD","fam","sex","Ancestry")
  keep_cols.adhd     <- c("interview_age_ADHD_min","ever_ADHD","ever_ADHD_str")
  keep_cols.comorbid <- c("IQ","mood","anxiety","externalizing","OCD","ODD","CD")
  keep_cols.prs      <- c("PRS",prs_measures)
  keep_cols.neurocog <- c("nihtbx_picvocab_fc_max","nihtbx_flanker_fc_max","nihtbx_list_fc_max","nihtbx_cardsort_fc_max","nihtbx_pattern_fc_max","nihtbx_picture_fc_max","nihtbx_reading_fc_max","nihtbx_fluidcomp_fc_max","nihtbx_cryst_fc_max","nihtbx_totalcomp_fc_max","pea_wiscv_tss_max")
  keep_cols.pc       <- c("PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
  keep_cols <- c(keep_cols.basic,keep_cols.adhd,keep_cols.comorbid,keep_cols.prs,keep_cols.neurocog,keep_cols.pc)
  sub_df <- sub_df[,keep_cols]

  #Merge demographic info with burden info
  analysis_df <- merge(x=sub_df,y=burden_df,by = "ID")
  rownames(analysis_df) <- analysis_df$ID
  #Set the current age to the last known interview date
  analysis_df$age <- analysis_df$interview_age_max

  #Return burden info
  return(analysis_df)
}

one_per_family_sample_targeted <- function(sub_df) {
  #Select one person per family, prioritizing ADHD females
  ids <- c()
  #loop through each family
  for (fam in unique(sub_df$fam)) {
    #Subset family based on sex and ADHD params
    fam_df <- sub_df[sub_df$fam == fam,]
    fam_df$sex_dx <- paste(fam_df$sex,fam_df$ever_ADHD,sep="_")
    #Extract a random family member that matches the criteria
    if (dim(fam_df)[1] == 1) {
      ids <- append(ids, fam_df[1,]$ID)
    } else if (dim(fam_df)[1] > 1) {
      for (sex_dx in c("F_ADHD","M_ADHD","M_non-ADHD","F_non-ADHD","other")) {
        if (!(sex_dx == "other")) {
          fam_df_sub <- fam_df[fam_df$sex_dx == sex_dx,]
          if (dim(fam_df_sub)[1] == 1) {
            ids <- append(ids, fam_df_sub[1,]$ID)
            break
          } else if (dim(fam_df_sub)[1] > 1) {
            ids <- append(ids, fam_df_sub[sample(rownames(fam_df_sub),size=1),]$ID)
            break
          }
        } else {
          ids <- append(ids, fam_df[sample(rownames(fam_df),size=1),]$ID)
          break
        }
      }
    } else {
      print(glue("ERROR: {fam} is empty"))
    }
  }
  sub_df <- sub_df[ids,]
  print(dim(sub_df))
  return(sub_df)
}

load_mri <- function(file,sub_df,type) {
  #Loads the MRI data
  mri <- read.csv(file)
  mri <- mri[mri$subjectkey %in% sub_df$ID,]
  if (type == "rsfmri") {
    #Exclude data with poor quality rsfmri
    mri <- mri[mri$imgincl_rsfmri_include == 1,]
    #List of highly correlated QC metrics for PCA
    pc_cols <- c("rsfmri_cor_ngd_scs_meanmn", "rsfmri_cor_ngd_scs_meanrot", "rsfmri_cor_ngd_scs_meantrans",
    "rsfmri_cor_ngd_scs_maxmn", "rsfmri_cor_ngd_scs_maxrot", "rsfmri_cor_ngd_scs_maxtrans")
    # Perform PCA
    res.pca <- prcomp(mri[,pc_cols],scale = TRUE)

    # Results for individuals
    res.ind <- get_pca_ind(res.pca)

    # Save PCs
    mri$rsfmri_pc1 <- res.ind$coord[,1]
    mri$rsfmri_pc2 <- res.ind$coord[,2]

    #Calculate quadratics of motion QC metrics
    mri$quad_rsfmri_cor_ngd_scs_stcgnvols <- mri$rsfmri_cor_ngd_scs_stcgnvols ^ 2
    mri$quad_rsfmri_pc1 <- mri$rsfmri_pc1 ^ 2
    mri$quad_rsfmri_pc2 <- mri$rsfmri_pc2 ^ 2
  } else if (type == "dmri") {
    #Exclude data with poor quality dmri
    mri <- mri[mri$imgincl_dmri_include == 1,]

    #Calculate quadratics for QC params
    mri$quad_dmri_dti_meanmotion <- mri$dmri_dti_meanmotion ^ 2
    mri$quad_dmri_dti_meantrans  <- mri$dmri_dti_meantrans  ^ 2
    mri$quad_dmri_dti_meanrot    <- mri$dmri_dti_meanrot    ^ 2
  } else {
    #Exclude data with poor quality smri
    mri <- mri[(mri$imgincl_t1w_include == 1) & (mri$imgincl_t2w_include == 1),]
  }
  mri$ID <- mri$subjectkey
  return(mri)
}

load_data <- function(input,pipeline,pipeline_filter,sex = "all",dx="ksads",opf=TRUE,exclude = TRUE,equal=FALSE) {
  #Determine sex suffix
  if (sex == "all") {
    sex_str <- ""
  } else {
    sex_str <- glue("_{sex}")
  }
  #Determine dx suffix
  if (dx == "ksads") {
    dx_str <- ""
  } else {
    dx_str <- glue("_{dx}")
  }
  #Determine file name
  file.name <- glue("results/{pipeline}/data_files/{pipeline_filter}{sex_str}{dx_str}.rds")

  #Check if the data file already exists
  if (file.exists(file.name)) {
    print(glue("Loading {file.name}..."))
    #Load file if it exists
    out_df <- readRDS(file.name)
  } else {
    #Generate the data file
    print(glue("Generating {file.name}..."))
    dir.create(glue("results/{pipeline}/data_files/"),showWarnings = FALSE, recursive = TRUE)
    #Load burdens
    out_df <- load_abcd_burden(input,pipeline,pipeline_filter,sex,dx=dx)
    #Select one sample per family
    if(opf == TRUE) {
      out_df <- one_per_family_sample_targeted(out_df)
    }

    #Save file
    saveRDS(out_df, file = file.name)
  }
  if (exclude) {
    #Limit cohort by removing subjects excluded for various reasons
    out_df <- out_df[!startsWith(out_df$ever_ADHD_str,"EXCLUDE"),]
  }
  if (equal) {
    #Equalize proportion of males and females in analysis groups (robustness check)
    equal_type = "reduce_males"
    file.name <- glue("results/{pipeline}/data_files/{pipeline_filter}{sex_str}{dx_str}_equal_{equal_type}.rds")
    if (file.exists(file.name)) {
      print(glue("Loading {file.name}..."))
      #Load file if it exists
      out_df <- readRDS(file.name)
    } else {
      # Limit to same number for all groups
      out_df$sex_ADHD <- paste(out_df$sex,out_df$ever_ADHD,sep="_")
      out_df <- out_df[out_df$ever_ADHD %in% c("ADHD","non-ADHD"),]
      if (equal_type == "reduce_males") {
        print("Reducing sample size of ADHD males to equalize proportions")
        #Add all females to output
        out_df2 <- out_df[out_df$sex_ADHD %in% c("F_ADHD","F_non-ADHD"),]
        #Add all non-ADHD males to output
        out_df2 <- rbind(out_df2,out_df[out_df$sex_ADHD == "M_non-ADHD",])
        # Calculate ADHD/non-ADHD ratio for females
        ratio <- sum(out_df$sex_ADHD == "F_ADHD")/sum(out_df$sex_ADHD == "F_non-ADHD")
        print(ratio)
        #Reduce the number of male ADHD by the given ratio
        out_df_tmp <- out_df[out_df$sex_ADHD == "M_ADHD",]
        out_df_tmp <- out_df_tmp[sample(rownames(out_df_tmp),size = round(ratio*nrow(out_df[out_df$sex_ADHD == "M_non-ADHD",]))),]
      } else if (equal_type == "reduce_females") {
        print("Reducing sample size of non-ADHD females to equalize proportions")
        #Add all ADHD cases to output
        out_df2 <- out_df[out_df$sex_ADHD %in% c("F_ADHD","M_ADHD"),]
        #Add all non-ADHD males to output
        out_df2 <- rbind(out_df2,out_df[out_df$sex_ADHD == "M_non-ADHD",])
        # Calculate non-ADHD/ADHD ratio for males
        ratio <- sum(out_df$sex_ADHD == "M_non-ADHD")/sum(out_df$sex_ADHD == "M_ADHD")
        #Reduce the number of female non-ADHD by the given ratio
        out_df_tmp <- out_df[out_df$sex_ADHD == "F_non-ADHD",]
        out_df_tmp <- out_df_tmp[sample(rownames(out_df_tmp),size = round(ratio*nrow(out_df[out_df$sex_ADHD == "F_ADHD",]))),]
      }

      out_df <- rbind(out_df2,out_df_tmp)
      saveRDS(out_df, file = file.name)
    }

  }
  #For case/control
  print("Final Data:")
  for (status in c("ADHD","non-ADHD")) {
    print(status)
    print(glue("\tN = {sum(out_df$ever_ADHD == status,na.rm = TRUE)}"))
    print(glue("\t%_Female = {sum(out_df$ever_ADHD == status & out_df$sex == 'F',na.rm = TRUE)/sum(out_df$ever_ADHD == status,na.rm = TRUE)*100}"))
    print(glue("\tAge (mean) = {mean(out_df[out_df$ever_ADHD == status,'age'],na.rm = TRUE)}"))
    print(glue("\tAge (std) = {sd(out_df[out_df$ever_ADHD == status,'age'],na.rm = TRUE)}"))
  }
  print("Total")
  print(glue("\tN = {nrow(out_df)}"))
  print(glue("\t%_Female = {sum(out_df$sex == 'F',na.rm = TRUE)/nrow(out_df)*100}"))
  print(glue("\tAge (mean) = {mean(out_df$age,na.rm = TRUE)}"))
  print(glue("\tAge (std) = {sd(out_df$age,na.rm = TRUE)}"))
  write.csv(out_df,"scripts/11_revisions/subject_df.csv")
  return(out_df)
}
