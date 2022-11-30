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
} else if (length(args) == 1) {
  pipeline = args[1]
  pipeline_filter = "bpd_psych"
  sex = "all"
} else if (length(args) == 2) {
  pipeline = args[1]
  pipeline_filter = args[2]
  sex = "all"
} else if (length(args) == 3) {
  pipeline = args[1]
  pipeline_filter = args[2]
  sex = args[3]
}

out_name <- "mediation/cognitive"


input <- "consensus"

#Load burden and subject info
sub_df <- load_data(input,pipeline,pipeline_filter,sex = sex)

if (sex != "all") {
  pipeline_filter <- paste(pipeline_filter,sex,sep="_")
}

print(pipeline_filter)
summarize_mediation <- function(results,param) {
  #Create a summary table for the mediation results
  if (is.null(results)) {
    #If no results, create an empty table to return
    smat <- matrix(nrow=1,ncol=4)
  } else {
    #Get original summary table
    res <- summary(results)
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

run_mediation <- function(df,mod_param,cov_params,burden) {
n_sim = 10000

#Scale mediator
df$cog <- scale(df[,mod_param],center=TRUE,scale=TRUE)[,]
#Remove NA values from cog metric
df <- df[!is.na(df[,mod_param]),]

#Define mediator model
model.M <- lm(as.formula(paste("cog ~",burden,"+", paste(cov_params,collapse = " + "),sep =  " ")), data = df)
#If model.M is not significant, then there is no justification for running mediation
#Define combined model
model.Y <- glm(as.formula(paste("ever_ADHD ~",burden,"+ cog +", paste(cov_params,collapse = " + "), sep =  " ")), data = df, family="binomial")
#Perform mediation
if (startsWith(burden, "Binary")) {
  results <- mediate(model.M, model.Y, treat=burden, mediator="cog",
                     treat.value = TRUE, control.value = FALSE,
                     boot=TRUE, sims=n_sim)
} else {
  results <- mediate(model.M, model.Y, treat=burden, mediator="cog",
                     boot=TRUE, sims=n_sim)
}
res <- summarize_mediation(results,mod_param)

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
res_sup[i,"param"] <- mod_param

#Indirect (M -> DV)
i = 2
ci <- confint(model.Y)
print(summary(model.Y))
res_sup[i,"Effect"] <- "Indirect (DV ~ M)"
# print(ci)
res_sup[i,"0.025_CI"]  <- ci["cog","2.5 %"]
res_sup[i,"0.975_CI"] <- ci["cog","97.5 %"]
res_sup[i,"Estimate"] <- coef(summary(model.Y))["cog","Estimate"]
res_sup[i,"p-value"] <- coef(summary(model.Y))["cog","Pr(>|z|)"]
res_sup[i,"param"] <- mod_param

res <- rbind(res, res_sup)
print(glue("Mediation model complete for param {mod_param}"))
res$burden <- burden
return(res)
}


#Set up parallel backend
myCluster <- makeCluster(28, # number of cores to use
                         type = "PSOCK") # type of cluster

registerDoParallel(myCluster)


#Define covariates to use in model: starting with the batch and SNP PCs AND site
cov_params <- c("IQ","BATCH_METHOD","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")

sub_df_opf <- sub_df

if (length(unique(sub_df_opf$sex)) == 2) {
  cov_params <- append("sex",cov_params)
}

burden_of_int <- c("Binary_200kb","Length","N_genes")

cog_list <- c("nihtbx_picvocab_fc","nihtbx_flanker_fc","nihtbx_list_fc","nihtbx_cardsort_fc","nihtbx_pattern_fc","nihtbx_picture_fc","nihtbx_reading_fc","IQ")

cnv_event <- "del"

print(glue("\ttesting {cnv_event} CNV types..."))
#Start timer for processing time
ptm <- proc.time()
outdir <- glue("results/{pipeline}/stats/{pipeline_filter}/{out_name}/{cnv_event}")
dir.create(outdir, recursive = TRUE,showWarnings = FALSE)

# Subset cog columns to only those of interest
print(glue("\trunning mediation analysis on {length(cog_list)} params..."))
# Loop through each cog param
results <- foreach(mod_param = rep(cog_list,length(burden_of_int)),burden = rep(burden_of_int,each = length(cog_list)),.combine = 'rbind',.export = ls(globalenv()),.packages = c("glue","stringr","mediation")) %dorng% {
  burden <- glue("{burden}_{cnv_event}")
  .GlobalEnv$burden <- burden
  .GlobalEnv$cov_params <- cov_params
  print(burden)
  #Run mediation analysis
  print(head(sub_df_opf))
  results <- run_mediation(sub_df_opf,mod_param,cov_params,burden)
  return(results)
}
#Write results
write.csv(results,glue("{outdir}/Cognitive_mediation_results.csv"))



start_flag = TRUE
for (cnv_event in cnv_types) {
    input <- glue("results/{pipeline}/stats/{pipeline_filter}/{out_name}/{cnv_event}/Cognitive_mediation_results.csv")
    if (file.exists(input)) {
      df <- read.csv(input)
      print(head(df))
      df$CNV <- cnv_event
      df$cog <- "cognition"
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
            df[df$param_burden == param,"Description"] <- "Significant"
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

write.csv(out_df,glue("results/{pipeline}/stats/{pipeline_filter}/{out_name}/COG_mediation_results.csv"))
