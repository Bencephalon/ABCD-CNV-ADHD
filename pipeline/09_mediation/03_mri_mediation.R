library(mediation)
library(glue)
library(tidyverse)
library(stringr)
library(foreach)
library(doParallel)
library(doRNG)
source("scripts/99_util/util_load_cnvs.R")

set.seed(2514)

#Load pipeline
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  pipeline = "consensus"
  pipeline_filter = "bpd_psych"
  sex = "all"
  specificity = "targeted"
} else if (length(args) == 1) {
  pipeline = args[1]
  pipeline_filter = "bpd_psych"
  sex = "all"
  specificity = "targeted"
} else if (length(args) == 2) {
  pipeline = args[1]
  pipeline_filter = args[2]
  sex = "all"
  specificity = "targeted"
} else if (length(args) == 3) {
  pipeline = args[1]
  pipeline_filter = args[2]
  sex = args[3]
  specificity = "targeted"
} else if (length(args) == 4) {
  pipeline = args[1]
  pipeline_filter = args[2]
  sex = args[3]
  specificity = "targeted"
} else if (length(args) == 5) {
  pipeline = args[1]
  pipeline_filter = args[2]
  sex = args[3]
  specificity = args[4]
}

input <- "consensus"

#Load burden and subject info
sub_df <- load_data(input,pipeline,pipeline_filter,sex = sex)

if (sex != "all") {
  pipeline_filter <- paste(pipeline_filter,sex,sep="_")
}

if (specificity == "general") {
  n = 10000
  #Load in params
  mri.desc <- read.csv("data/mri/processed/sMRI.desc.csv",row.names=1)
  #Grepl params
  params_of_int <- mri.desc[grepl("smri_thick_cdk_",mri.desc$col) | grepl("smri_area_cdk_",mri.desc$col) | grepl("smri_vol_cdk_",mri.desc$col),"col"]
} else {
  n = 10000
  params_of_int <- c("smri_thick_cdk_meanlh","smri_thick_cdk_meanrh","smri_thick_cdk_mean",
                     "smri_area_cdk_totallh","smri_area_cdk_totalrh","smri_area_cdk_total",
                     "smri_vol_cdk_totallh" ,"smri_vol_cdk_totalrh" ,"smri_vol_cdk_total",
                     "smri_vol_scs_amygdalalh","smri_vol_scs_amygdalarh",
                     "smri_vol_scs_caudatelh","smri_vol_scs_caudaterh",
                     "smri_vol_scs_putamenlh","smri_vol_scs_putamenrh",
                     "smri_vol_scs_hpuslh","smri_vol_scs_hpusrh")
}

out_name <- glue("mediation_{specificity}")


n_sim = n

summarize_mediation <- function(results,param) {
  #Create a summary table for the mediation results
  if (is.null(results)) {
    #If no results, create an empty table to return
    smat <- matrix(nrow=1,ncol=4)
  } else {
    #Get original summary table
    res <- summary(results)
    #https://stackoverflow.com/questions/41582486/how-to-convert-r-mediation-summary-to-data-frame
    #Check model type
    isLinear.y <- ((class(res$model.y)[1] %in% c("lm", "rq")) ||
                     (inherits(res$model.y, "glm") && res$model.y$family$family ==
                        "gaussian" && res$model.y$family$link == "identity") ||
                     (inherits(res$model.y, "survreg") && res$model.y$dist ==
                        "gaussian"))

    printone <- !res$INT && isLinear.y
    #Compile results
    if (printone) {
      smat <- c(res$d1, res$d1.ci, res$d1.p)
      smat <- rbind(smat, c(res$z0, res$z0.ci, res$z0.p))
      smat <- rbind(smat, c(res$tau.coef, res$tau.ci, res$tau.p))
      smat <- rbind(smat, c(res$n0, res$n0.ci, res$n0.p))
    } else {
      smat <- c(res$d0, res$d0.ci, res$d0.p)
      smat <- rbind(smat, c(res$d1, res$d1.ci, res$d1.p))
      smat <- rbind(smat, c(res$z0, res$z0.ci, res$z0.p))
      smat <- rbind(smat, c(res$z1, res$z1.ci, res$z1.p))
      smat <- rbind(smat, c(res$tau.coef, res$tau.ci, res$tau.p))
      smat <- rbind(smat, c(res$n0, res$n0.ci, res$n0.p))
      smat <- rbind(smat, c(res$n1, res$n1.ci, res$n1.p))
      smat <- rbind(smat, c(res$d.avg, res$d.avg.ci, res$d.avg.p))
      smat <- rbind(smat, c(res$z.avg, res$z.avg.ci, res$z.avg.p))
      smat <- rbind(smat, c(res$n.avg, res$n.avg.ci, res$n.avg.p))
    }
  }
  #Define data columns
  colnames(smat) <- c("Estimate", "0.025_CI", "0.975_CI", "p-value")
  #Convert to data frame
  smat <- as.data.frame(smat)
  #Rename rows
  rownames(smat) <- 1:dim(smat)[1]
  #Define what results each row is displaying
  if (is.null(results)) {
    smat$Effect <- "Invalid"
  } else if (printone) {
    smat$Effect <- c("ACME", "ADE", "Total Effect", "Prop. Mediated")
  } else {
    smat$Effect <- c("ACME (control)", "ACME (treated)",
                        "ADE (control)", "ADE (treated)", "Total Effect",
                        "Prop. Mediated (control)", "Prop. Mediated (treated)",
                        "ACME (average)", "ADE (average)", "Prop. Mediated (average)")
  }
  #Define the param that was tested to keep track after rows combined
  smat$param <- param
  #Return the summary table
  return(smat)
}

run_mediation <- function(df,mri_param,cov_params,burden) {
  #Scale mediator
  df$MRI <- scale(df[,mri_param],center=TRUE,scale=TRUE)[,]
  #Remove NA values from MRI metric
  df <- df[!is.na(df[,mri_param]),]

  #Define mediator model
  model.M <- lm(as.formula(paste("MRI ~",burden,"+", paste(cov_params,collapse = " + "),sep =  " ")), data = df)
  #Define combined model
  model.Y <- glm(as.formula(paste("ever_ADHD ~",burden,"+ MRI +", paste(cov_params,collapse = " + "), sep =  " ")), data = df, family="binomial")
  # model.Y <- glmer(as.formula(paste("ever_ADHD ~",burden,"+ MRI +", paste(cov_params,collapse = " + "), sep =  " ")), data = df, family="binomial")
  #Perform mediation
  if (startsWith(burden, "Binary")) {
    results <- mediate(model.M, model.Y, treat=burden, mediator="MRI",
                       treat.value = TRUE, control.value = FALSE,
                       boot=TRUE, sims=n_sim)
  } else {
    results <- mediate(model.M, model.Y, treat=burden, mediator="MRI",
                       boot=TRUE, sims=n_sim)
  }

  res <- summarize_mediation(results,mri_param)

  res_sup <- data.frame(matrix(nrow=2,ncol=length(colnames(res))))
  colnames(res_sup) <- colnames(res)

  if (startsWith(burden, "Binary")) {
     b_match = glue("{burden}TRUE")
     p_match = "Pr(>|t|)"
  } else {
     b_match = burden
     p_match = "Pr(>|z|)"
  }
  #Indirect (IV -> M)
  i = 1
  ci <- confint(model.M)
  res_sup[i,"Effect"] <- "Indirect (M ~ IV)"
  res_sup[i,"0.025_CI"]  <- ci[b_match,"2.5 %"]
  res_sup[i,"0.975_CI"] <- ci[b_match,"97.5 %"]
  res_sup[i,"Estimate"] <- coef(summary(model.M))[b_match,"Estimate"]
  res_sup[i,"p-value"] <- coef(summary(model.M))[b_match,"Pr(>|t|)"]
  res_sup[i,"param"] <- mri_param

  #Indirect (M -> DV)
  i = 2
  ci <- confint(model.Y)
  res_sup[i,"Effect"] <- "Indirect (DV ~ M)"
  # print(ci)
  res_sup[i,"0.025_CI"]  <- ci["MRI","2.5 %"]
  res_sup[i,"0.975_CI"] <- ci["MRI","97.5 %"]
  res_sup[i,"Estimate"] <- coef(summary(model.Y))["MRI","Estimate"]
  res_sup[i,"p-value"] <- coef(summary(model.Y))["MRI","Pr(>|z|)"]
  res_sup[i,"param"] <- mri_param

  res <- rbind(res, res_sup)
  res$burden <- burden
  print(res)
  return(res)
}


#Set up parallel backend
myCluster <- makeCluster(28, # number of cores to use
                         type = "PSOCK") # type of cluster

registerDoParallel(myCluster)

################################# Load sMRI #################################
smri <- load_mri("data/mri/processed/sMRI.csv",sub_df,"smri")
# smri_cols <- colnames(smri)[str_detect(colnames(smri),"smri_.*_cdk_.*")]
smri_cols <- colnames(smri)

#Define covariates to use in model: starting with the batch and SNP PCs AND site
cov_params <- c("IQ","interview_age","BATCH_METHOD","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

#Select all subjects with sMRI
sub_df <- sub_df[sub_df$ID %in% smri$ID,]

smri_df <- merge(sub_df,smri,by = c("ID","sex"))

#Normalize all metrics by relevant metric
normalized_cols <- params_of_int[params_of_int %in% smri_cols]

#Noramlize cortical volume
normalized_cols_sub <- normalized_cols[grepl("smri_vol_scs",normalized_cols)]
normalized_cols_sub <- normalized_cols_sub[normalized_cols_sub != "smri_vol_cdk_total"]
smri_df[,normalized_cols_sub] <- smri_df[,normalized_cols_sub]/smri_df$smri_vol_cdk_total

if (length(unique(sub_df$sex)) == 2) {
  cov_params <- append("sex",cov_params)
}


burden_of_int <- c("Binary_200kb","Length","N_genes")

cnv_event <- "del"
#Loop through each MRI type
mri_type <- "sMRI"
mri_cols <- smri_cols
df <- smri_df
cov_params_model <- cov_params
print(dim(mri_cols))
mri_cols <- mri_cols[mri_cols %in% colnames(df[,!(colSums(is.na(df))>0)])]
mri_cols <- mri_cols[apply(df[,mri_cols], 2, function(x) length(unique(x))) > 1]
mri_cols <- mri_cols[mri_cols %in% params_of_int]

for (cnv_event in cnv_types) {
  print(glue("\ttesting {cnv_event} CNV types..."))
  #Start timer for processing time
  ptm <- proc.time()
  outdir <- glue("results/{pipeline}/stats/{pipeline_filter}/{out_name}/{cnv_event}")
  dir.create(outdir, recursive = TRUE,showWarnings = FALSE)

  # Subset MRI columns to only those of interest
  print(glue("\trunning mediation analysis on {length(mri_cols)} params..."))
  # Loop through each MRI param
  results <- foreach(mri_param = rep(mri_cols,length(burden_of_int)),burden = rep(burden_of_int,each = length(mri_cols)),.combine = 'rbind',.export = ls(globalenv()),.packages = c("glue","stringr","mediation")) %dorng% {
    burden <- glue("{burden}_{cnv_event}")
    .GlobalEnv$burden <- burden
    .GlobalEnv$cov_params <- cov_params
    #Run moderated-mediation analysis
    results <- run_mediation(df,mri_param,cov_params,burden)
    return(results)
  }
  #Write results
  write.csv(results,glue("{outdir}/{mri_type}_mediation_results.csv"))
  print(proc.time() - ptm)
  print("-----")
}

#Load in description files
mri_params.desc <- read.csv("data/mri/processed/DTI.desc.csv")
mri_params.desc <- rbind(mri_params.desc, read.csv("data/mri/processed/sMRI.desc.csv"))
mri_params.desc <- rbind(mri_params.desc, read.csv("data/mri/processed/rsfMRI.desc.csv"))

start_flag = TRUE
for (cnv_event in cnv_types) {
  for (mri_event in c("sMRI","DTI","rs-fMRI")) {
    input <- glue("results/{pipeline}/stats/{pipeline_filter}/{out_name}/{cnv_event}/{mri_event}_mediation_results.csv")
    if (file.exists(input)) {
      # df <- read_csv(input)
      df <- read.csv(input)
      print(head(df))
      df$CNV <- cnv_event
      df$MRI <- mri_event
      df$p_ACME = -1
      df$p_Total = -1
      df$param_burden <- paste(df$param,df$burden,sep="-")
      df$Description <- ""

      for (param in unique(df$param_burden)) {
        p_ACME <- df$p.value[df$param_burden == param & df$Effect %in% c("ACME","ACME (average)","Invalid")]
        p_Total <- df$p.value[df$param_burden == param & df$Effect %in% c("Total Effect")]
        df[df$param_burden == param,"p_ACME"] <- p_ACME
        df[df$param_burden == param,"p_Total"] <- p_Total
        if (!is.na(p_ACME)) {
          if (p_ACME <= 0.05 & p_Total <= 0.05) {
            print(str_split(param,"-")[[1]][1])
            print(mri_params.desc$description[mri_params.desc$col == str_split(param,"-")[[1]]][1])
            df[df$param_burden == param,"Description"] <- mri_params.desc$description[mri_params.desc$col == str_split(param,"-")[[1]][1]]
          }
        }
      }


      if (start_flag) {
        out_df <- df
        start_flag = FALSE
      } else {
        print("Binding outputs...")
        out_df <- rbind(out_df,df)
      }
    } else {
      print(glue("{input} does not exist"))
    }

  }
}

write.csv(out_df,glue("results/{pipeline}/stats/{pipeline_filter}/{out_name}/MRI_mediation_results.csv"))
