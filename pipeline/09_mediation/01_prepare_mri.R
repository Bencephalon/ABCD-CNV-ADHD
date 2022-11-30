library(glue)
library(readxl)

setwd("/data/NCR_SBRB/jungbt/CNV/ABCD")

convert_to_integer <- function(col) {
  # print(head(col))
  tryCatch(
    {
      return(as.numeric(col))
    },
    warning = function(w)
    {
      return(col)
    }
  )
}

load_abcd <- function(loc) {
  df      <- read.table(loc,sep="\t",header=TRUE)
  df.desc <- data.frame(col = colnames(df),description = as.character(df[1,]))
  df      <- df[2:dim(df)[1],]
  for (col in colnames(df)) {
    df[,col] <- convert_to_integer(df[,col])
  }
  df$interview_age <- df$interview_age/12
  df$unique_id <- paste(df$subjectkey,df$eventname,sep="_")

  return(list(df,df.desc))
}

load_mri_data <- function(loc) {
  out <- load_abcd(loc)
  out[[1]] <- out[[1]][out[[1]]$unique_id %in% qc$unique_id,!(colnames(out[[1]]) %in% c("study_cohort_name","collection_title","eventname","src_subject_id"))]
  return(out)
}


merge_mri <- function(loc.1,loc.2,loc.3 = NULL, out_pref = "",site = site_df) {
  #Open mri.1 file
  out <- load_mri_data(loc.1)
  mri.1 <- out[[1]]
  mri.1.desc <- out[[2]]

  #Open mri.2 file
  out <- load_mri_data(loc.2)
  mri.2 <- out[[1]]
  mri.2.desc <- out[[2]]

  #Combine MRI files
  mri <- merge(mri.1,mri.2,by = c("unique_id","subjectkey","interview_age","sex"))
  #Combined MRI file descriptions
  mri.desc <- rbind(mri.1.desc,mri.2.desc)


  if (!is.null(loc.3)) {
    #Open mri.2 file
    out <- load_mri_data(loc.3)
    mri.3 <- out[[1]]
    mri.3.desc <- out[[2]]

    #Combine MRI files
    mri <- merge(mri,mri.3,by = c("unique_id","subjectkey","interview_age","sex"))
    #Combined MRI file descriptions
    mri.desc <- rbind(mri.desc,mri.3.desc)
  }

  mri <- merge(qc,mri,by = c("unique_id","subjectkey","interview_age","sex"))

  #Add site
  mri <- merge(mri,site,by = "subjectkey")

  #Write files
  write.csv(mri,glue("data/mri/processed/{out_pref}.csv"))
  write.csv(mri.desc,glue("data/mri/processed/{out_pref}.desc.csv"))

}

#Open QC file
qc.loc  <- "data/mri/abcd_imgincl01.txt"
out <- load_abcd(qc.loc)
qc <- out[[1]]
qc.desc <- out[[2]]
#Remove unnecessary columns
exclude <- c("collection_title","study_cohort_name","interview_date","src_subject_id","collection_id","abcd_imgincl01_id","dataset_id","visit")
qc      <- qc[, !names(qc) %in% exclude]
#Define QC columns
qc_inc  <- c("imgincl_t1w_include","imgincl_t2w_include")
#Calculate entries where all MRI data is usable
qc$all_include <- apply(qc[,qc_inc], 1, FUN = min)
#Calculate entries where no MRI data is usable
qc$all_exclude <- 1 - apply(qc[,qc_inc], 1, FUN = max)
#Remove entries where no MRI data is usable
qc <- qc[qc$all_include == 1,]

ages <- c(8, 8.5 ,9, 9.5, 10, 10.5, 11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5)
counts <- c()
for (i in ages) {
  out <- length(qc$interview_age[qc$interview_age >= i-0.5 & qc$interview_age < i + 0.5])
  print(glue("{i} years: {out}"))
  counts <- append(counts,out)
}

age.sel <- 10
print(age.sel)

qc$interview_age_net <- abs(age.sel - qc$interview_age)
qc <- qc[order(qc$interview_age_net),]
#Keep only the assessment closest to the target age
qc <- qc[!duplicated(qc$subjectkey),]
# hist(qc$interview_age_net)
print(glue("Mean Age Offset: {mean(qc$interview_age_net)}"))
print(glue("SD Age Offset: {sd(qc$interview_age_net)}"))

mean_age <- mean(qc$interview_age)
sd_age   <- sd(qc$interview_age)
print(glue("Mean Age: {mean_age}"))
print(glue("SD Age: {sd_age}"))

#Remove assessments outside of 3 SD
qc <- qc[qc$interview_age >= (mean_age - 3*sd_age) & qc$interview_age <= (mean_age + 3*sd_age),]

# Load in subject information
sub_df.loc <- file.path("results","consensus","phenotypes.csv")
sub_df <- read.csv(sub_df.loc)

#Limit QC data frame to subject data frame
qc <- qc[qc$subjectkey %in% sub_df$ID,]

print(glue("CNV Subjects Lacking MRI: {length(setdiff(sub_df$ID,qc$subjectkey))}"))

for (crit in qc_inc) {
  print(glue("\tSubjects w/ {crit}: {sum(qc[,crit])}"))
}

#Write to output
dir.create("data/mri/processed/", recursive = TRUE,showWarnings = FALSE)

#Load site data
site_df <- read_excel("data/mri/ABCD_metainfo.xlsx")
site_df <- site_df[!duplicated(site_df$ID...1),]
site_df <- site_df[,c("alt_ID","Scan_Location")]
colnames(site_df) <- c("subjectkey","site")

#Combine sMRI metrics
merge_mri("data/mri/abcd_smrip10201.txt","data/mri/abcd_smrip20201.txt","data/mri/abcd_smrip30201.txt",out_pref = "sMRI")

#Combine rs-fMRI metrics
merge_mri("data/mri/abcd_betnet02.txt","data/mri/mrirscor02.txt",out_pref = "rsfMRI")

#Combine DTI metrics
merge_mri("data/mri/abcd_dti_p101.txt","data/mri/abcd_dti_p201.txt",out_pref = "DTI")
