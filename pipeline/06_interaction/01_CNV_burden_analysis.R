# Libaries
library(tidyverse)
library(glue)
library(stringr)
library(emmeans)
source("scripts/99_util/util_load_cnvs.R")

pdf.options(encoding='ISOLatin2.enc')
set.seed(2514)

#Load pipeline
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  pipeline = "consensus"
  pipeline_filter = "bpd_psych"
  measure = "CNV"
  disorder = "ever_ADHD"
  dx = "ksads"
} else if (length(args) == 1) {
  pipeline = args[1]
  pipeline_filter = "bpd_psych"
  measure = "CNV"
  disorder = "ever_ADHD"
  dx = "ksads"
} else if (length(args) == 2) {
  pipeline = args[1]
  pipeline_filter = args[2]
  measure = "CNV"
  disorder = "ever_ADHD"
  dx = "ksads"
} else if (length(args) == 3) {
  pipeline = args[1]
  pipeline_filter = args[2]
  measure = args[3]
  disorder = "ever_ADHD"
  dx = "ksads"
} else if (length(args) == 4) {
  pipeline = args[1]
  pipeline_filter = args[2]
  measure = args[3]
  disorder = args[4]
  dx = "ksads"
} else if (length(args) == 5) {
  pipeline = args[1]
  pipeline_filter = args[2]
  measure = args[3]
  disorder = args[4]
  dx = args[5]
}


#For female-only model assign flags that allow for removal of sex terms
sex_flag <- FALSE
sex_t <- " + sex"


input <- "consensus"

#Load the data using the util helper functions
if (disorder == "ever_ADHD") {
  analysis_df <- load_data(input,pipeline,pipeline_filter,dx = dx,exclude = TRUE)
} else {
  analysis_df <- load_data(input,pipeline,pipeline_filter,dx = dx,exclude = FALSE)
}

#Rename the pipeline based on modifiers
if (!(dx == "ksads")) {
  pipeline_filter <- glue("{pipeline_filter}_{dx}")
}

#Define output columns
est_cols <- c("Estimate", "Std. Error","z value","P","CI_0.025","CI_0.975","burden")

#Define burden/PRS measures for analysis
if (measure == "CNV") {
  burden_types <- c("Binary_50kb","Binary_100kb","Binary_200kb","N_CNVs","N_genes","Length","LOEUF","pHI","pTS")
} else if (measure == "PRS") {
  analysis_df <- analysis_df[!is.na(analysis_df$PRS),]
  burden_types <- c("Pt_5e.08","Pt_1e.05","Pt_0.0001","Pt_0.001","Pt_0.01","Pt_0.05","Pt_0.1","Pt_0.2","Pt_0.5")
}



get_ors <- function(model,burden,m_type) {
  if (burden %in% c("Binary_50kb","Binary_100kb","Binary_200kb","Binary_Patho")){
    #Run simple emmeans test
    em <- emmeans(model,specs = c("sex",m_type))
    posthoc <- pairs(em,reverse= TRUE, simple = "each")
    out <- as.data.frame(posthoc[[2]])
    #Calculate OR
    out$OR <- exp(out$estimate)
    #Calculate OR CI
    out[,c("CI_lo","CI_hi")] <- exp(as.data.frame(confint(posthoc)[[2]])[,c("asymp.LCL","asymp.UCL")])
    # Format outputs
    out <- out[,!(colnames(out) %in% c("contrast","SE","df","z.ratio","estimate"))]
  } else {
    #Run simple emtrends test
    em <- emtrends(model,specs = c("sex",m_type),var = m_type)
    out <- as.data.frame(em)
    #Calculate OR
    out$OR <- exp(out[,glue("{m_type}.trend")])
    #Calculate OR CI
    out[,c("CI_lo","CI_hi")] <- exp(out[,c("asymp.LCL","asymp.UCL")])
    out$p.value <- test(em)$p.value
  }
  #Assign sex
  out$sex <- c("female","male")
  #Assign burden measure
  out$burden <- burden
  return(out[c("sex","OR","CI_lo","CI_hi","p.value","burden")])
}

run_model <- function(df, f, out_pref, out_dir,sex_t = " + sex",col_int = NULL,m_type = "") {
  dir.create(file.path(out_dir,"csvs"), recursive = TRUE,showWarnings = FALSE)
  df <- data.frame(df)

  models <- c()
  #Construct formula
  constant_covariates <- " + IQ + BATCH_METHOD + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
  f <- paste(f,constant_covariates,sep="")

  #Store results
  results_table <- data.frame(matrix(nrow = 0,ncol=length(est_cols)))
  colnames(results_table) <- est_cols

  #Store ORs
  or_cols <- c("sex","OR","CI_lo","CI_hi","p.value","burden")
  or_table <- data.frame(matrix(nrow = 0,ncol=length(or_cols)))
  colnames(or_table) <- or_cols

  #Loop through each burden metric
  for (burden in burden_types) {
      #Define type of model to run
      family <- "binomial"
      #Redefine specific burden metric to generic del and dup variable
      if (measure == "CNV") {
        df[,m_type] <- df[,glue("{burden}_{m_type}")]
      } else if (measure == "PRS") {
        df[,m_type] <- df[,glue("{burden}")]
      }
      #Run model
      formula <- as.formula(glue(f))
      model <- glm(formula, data = df, family = family)

      #Generate odds ratios
      results <- data.frame(cbind(coef(summary(model)),confint.default(model)))
      results$burden <- burden
      colnames(results) <- est_cols
      results_table <- rbind(results_table, results)

      #Get ORs
      if (col_int == "int") {
        or_table <- rbind(or_table,get_ors(model,burden,m_type))
      }

      print("===================================================")
  }
  #Save model results
  if (!is.null(col_int)) {
    if (col_int %in% c("del","dup","both","PRS")) {
      write.csv(results_table, file.path(out_dir,"csvs",paste(out_pref,"_est.csv",sep="")))
    } else if (col_int == "int") {
      write.csv(results_table, file.path(out_dir,"csvs",paste(out_pref,"_sexint_est.csv",sep="")))
      write.csv(or_table     , file.path(out_dir,"csvs",paste(out_pref,"_OR.csv",sep="")))
    }
  }
}




for (rm_low_iq in c(TRUE,FALSE)) {
  #Check whether to remove the lower IQ individuals
  if (rm_low_iq) {
    out_name <- "burden_IQ"
    analysis_df_tmp <- analysis_df[analysis_df$IQ_unscaled >= 3 & !is.na(analysis_df$IQ_unscaled),]
  } else {
    out_name <- "burden"
    analysis_df_tmp <- analysis_df
  }
  #Create output dir
  outdir <- file.path("results",pipeline,"stats",pipeline_filter,out_name,disorder,measure)
  print(outdir)
  dir.create(file.path(outdir,"csvs"), recursive = TRUE,showWarnings = FALSE)

  #Define regression measure
  if (measure == "CNV") {
    measures <- c("del","dup")
  } else if (measure == "PRS") {
    measures <- c("PRS")
  }
  #Loop through each regression measure
  for (m in measures) {
      print(m)
      #Redefine pHI/pTS as necessary
      if (m == "del") {
        analysis_df_tmp$pHI_pTS_del <- analysis_df_tmp$pHI_del
      } else if (m == "dup") {
        analysis_df_tmp$pHI_pTS_dup <- analysis_df_tmp$pTS_dup
      } else {
        analysis_df_tmp$pHI_pTS_both <- 0
      }
      #Run model without sex interaction
      print("Main")
      out <- glue("{m}")
      f <- sprintf("{disorder} ~ %s{sex_t} ",m)
      models   <- run_model(analysis_df_tmp, f, out, outdir, sex_t = sex_t, col_int = m,m_type = m)
      #Run model with sex interaction
      print("Interaction")
      out <- glue("{m}_sexI")
      f <- sprintf("{disorder} ~ %s{sex_t} + %s*sex ",m,m)
      models   <- run_model(analysis_df_tmp, f, out, outdir, sex_t = sex_t, col_int = "int",m_type = m)
  }
  #For main and interaction
  for (effect in c("main","interaction")) {
    #Make out summary df
    summary_df <- data.frame(matrix(nrow=length(burden_types)*length(measures),ncol=length(est_cols) + 2)) #nrow =
    colnames(summary_df) <- c("measure","burden_type",est_cols)
    summary_df$measure <- rep(measures,each=length(burden_types))
    summary_df$burden_type <- rep(burden_types,length(measures))
    #Loop through each model output and extract the necessary results
    for (row in 1:nrow(summary_df)) {
      #Load file
      if (effect == "main") {
        df <- read.csv(glue("{outdir}/csvs/{summary_df[row,'measure']}_est.csv"),row.names = 1)
        summary_df[row,est_cols] <- df[df$burden == summary_df[row,"burden_type"] & startsWith(rownames(df),summary_df[row,"measure"]),]
      } else {
        df <- read.csv(glue("{outdir}/csvs/{summary_df[row,'measure']}_sexI_sexint_est.csv"),row.names = 1)
        summary_df[row,est_cols] <- df[df$burden == summary_df[row,"burden_type"] & startsWith(rownames(df),summary_df[row,"measure"]) & grepl("sexM",rownames(df)),]
      }
    }
    # Remove dup pHI and del pTS
    summary_df <- summary_df[!(summary_df$burden == "pHI" & summary_df$measure == "dup"),]
    summary_df <- summary_df[!(summary_df$burden == "pTS" & summary_df$measure == "del"),]

    #Apply FDR correction
    summary_df$P_adjust <- p.adjust(summary_df$P,method = "fdr")
    #Save results
    summary_df <- summary_df[c("measure","Estimate","Std. Error","CI_0.025","CI_0.975","z value","P","P_adjust","burden")]
    write.csv(summary_df,glue("{outdir}/summary_{effect}.csv"))
  }
}
