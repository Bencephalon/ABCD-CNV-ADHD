# Libaries
library(tidyverse)
library(glue)
library(stringr)

source("scripts/99_util/util_load_cnvs.R")

pdf.options(encoding='ISOLatin2.enc')
set.seed(2514)

#Load pipeline
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  pipeline = "consensus"
  pipeline_filter = "bpd_psych"
  dx = "ksads"
} else if (length(args) == 1) {
  pipeline = args[1]
  pipeline_filter = "bpd_psych"
  dx = "ksads"
} else if (length(args) == 2) {
  pipeline = args[1]
  pipeline_filter = args[2]
  dx = "ksads"
} else if (length(args) == 3) {
  pipeline = args[1]
  pipeline_filter = args[2]
  dx = args[3]
}


#For female-only model assign flags that allow for removal of sex terms
sex_flag <- FALSE
sex_t <- " + sex"

input <- "consensus"

#Load in data
analysis_df <- load_data(input,pipeline,pipeline_filter,dx = dx)

if (!(dx == "ksads")) {
  pipeline_filter <- glue("{pipeline_filter}_{dx}")
}

est_cols <- c("Estimate", "Std. Error","z value","P","CI_0.025","CI_0.975","burden")
#Define moderators to test
moderators <- c("PRS", "nihtbx_picvocab_fc", "nihtbx_flanker_fc", "nihtbx_list_fc", "nihtbx_cardsort_fc", "nihtbx_pattern_fc", "nihtbx_picture_fc", "nihtbx_reading_fc", "nihtbx_fluidcomp_fc", "nihtbx_cryst_fc", "nihtbx_totalcomp_fc", "strp_scr_acc_all", "IQ")
#Define CNV measures to test
burden_types <- c("Binary_200kb","Length","N_genes")




run_model <- function(df, f, out_pref, out_dir,sex_t = " + sex",col_int = NULL,m_type = "") {
  dir.create(file.path(out_dir,"csvs"), recursive = TRUE,showWarnings = FALSE)
  df <- data.frame(df)

  models <- c()
  #Construct formula
  constant_covariates <- " + IQ + BATCH_METHOD + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
  f <- paste(f,constant_covariates,sep="")

  #Define burden metrics
  results_table <- data.frame(matrix(nrow = 0,ncol=length(est_cols)))
  colnames(results_table) <- est_cols

  #Loop through each burden metric
  for (burden in burden_types) {
      #Define type of model to run
      family <- "binomial"
      #Redefine specific burden metric to generic del and dup variable
      df[,m_type] <- df[,glue("{burden}_{m_type}")]
      #Run model
      formula <- as.formula(glue(f))
      model <- glm(formula, data = df, family = family)
      #Generate odds ratios
      results <- data.frame(cbind(coef(summary(model)),confint.default(model)))
      results$burden <- burden
      colnames(results) <- est_cols
      results_table <- rbind(results_table, results)
      print("===================================================")
  }
  #Save model results
  if (!is.null(col_int)) {
    if (col_int %in% c("del","dup","both","PRS")) {
      write.csv(results_table,file.path(out_dir,"csvs",paste(out_pref,"_est.csv",sep="")))
    } else if (col_int == "int") {
      write.csv(results_table,file.path(out_dir,"csvs",paste(out_pref,"_sexint_est.csv",sep="")))
    }
  }
}


out_name <- "moderation"
analysis_df_tmp <- analysis_df
measures <- c("del")

#Loop through each CNV measure
for (measure in measures) {
  #Loop through each potential moderator
  for (mod in moderators) {
    #Create output dirs
    outdir <- file.path("results",pipeline,"stats",pipeline_filter,out_name,mod)
    print(outdir)
    dir.create(file.path(outdir,"csvs"), recursive = TRUE,showWarnings = FALSE)
    #Run models with various interactions
    print("Model 1")
    #Model 1: Maximum Model: All interactions
    out <- glue("{measure}_model1")
    f <- sprintf(glue("ever_ADHD ~ {measure}*{mod}*sex + {measure}*{mod} + {measure}*sex + {mod}*sex "))
    print(f)
    models   <- run_model(analysis_df_tmp, f, out, outdir, col_int = "int",m_type = measure)
    print("Model 2")
    #Model 2: All 2-way interactions only
    out <- glue("{measure}_model2")
    f <- sprintf(glue("ever_ADHD ~ {measure}*{mod} + {measure}*sex + {mod}*sex "))
    models   <- run_model(analysis_df_tmp, f, out, outdir, col_int = "int",m_type = measure)
    print("Model 3")
    #Model 3: measure*sex and moderator*sex only
    out <- glue("{measure}_model3")
    f <- sprintf(glue("ever_ADHD ~ {measure}*sex + {mod}*sex "))
    models   <- run_model(analysis_df_tmp, f, out, outdir, col_int = "int",m_type = measure)
    print("Model 4")
    #Model 4: measure*moderator and moderator*sex only
    out <- glue("{measure}_model4")
    f <- sprintf(glue("ever_ADHD ~ {measure}*{mod} + {measure}*sex "))
    models   <- run_model(analysis_df_tmp, f, out, outdir, col_int = "int",m_type = measure)
    print("Model 5")
    #Model 5: measure*sex and moderator cov only
    out <- glue("{measure}_model5")
    f <- sprintf(glue("ever_ADHD ~ {measure}*sex + {mod} "))
    models   <- run_model(analysis_df_tmp, f, out, outdir, col_int = "int",m_type = measure)
    print("Model 6")
    #Model 6: moderator*sex and measure cov only
    out <- glue("{measure}_model6")
    f <- sprintf(glue("ever_ADHD ~ {mod}*sex + {measure} "))
    models   <- run_model(analysis_df_tmp, f, out, outdir, col_int = "int",m_type = measure)
  }
  #Make summary table
  #Load in all relevant files
  print(glue("{outdir}/csvs/*.csv"))
  files <- Sys.glob(glue("results/{pipeline}/stats/{pipeline_filter}/{out_name}/*/csvs/*.csv"))
  print(glue("Files: {length(files)}"))

  #Create dataframe to summarize results
  col_names <- c("burden","measure","model","n_interactions","n_sig",glue("p_del:sexM:mod"),glue("mod:sexM"),glue("p_del:mod"),glue("p_del:sexM"))
  summary_df <- data.frame(matrix(nrow=0,ncol=length(col_names)))
  colnames(summary_df) <- col_names

  #Loop through each result
  for (file in files) {
    #Read in the results file
    df <- read.csv(file)
    measure <- basename(dirname(dirname(file)))
    model <- str_split(basename(file),"_")[[1]][2]
    out_df <- data.frame(matrix(nrow=length(burden_types),ncol=length(col_names)))
    #Format results file
    rownames(out_df) <- burden_types
    colnames(out_df) <- col_names
    out_df$measure <- measure
    out_df$burden <- burden_types
    out_df$model <- model

    for (burden in burden_types) {
      #Search rows for ":" (n_interactions)
      out_df[burden,"n_interactions"] <- nrow(df[df$burden == burden & grepl(":",df$X),])
      #Determine sig of all relevant rows (n_sig)
      out_df[burden,"n_sig"] <- nrow(df[df$burden == burden & grepl(":",df$X) & df$P <= 0.05,])
      if (nrow(out_df) > 0) {
        #Define measure:sex:CNV
        tmp_df <- df[grepl(measure,df$X) & grepl("sexM",df$X) & grepl("del",df$X) & df$burden == burden,]
        if (nrow(tmp_df) == 1) {
          out_df[burden,glue("p_del:sexM:mod")] <- tmp_df$P
        }
        #Define measure:sex
        tmp_df <- df[grepl(measure,df$X) & grepl("sexM",df$X) & !grepl("del",df$X) & df$burden == burden,]
        if (nrow(tmp_df) == 1) {
          out_df[burden,glue("mod:sexM")] <- tmp_df$P
        }
        #Define measure:CNV
        tmp_df <- df[grepl(measure,df$X) & !grepl("sexM",df$X) & grepl("del",df$X) & df$burden == burden,]
        if (nrow(tmp_df) == 1) {
          out_df[burden,glue("p_del:mod")] <- tmp_df$P
        }
        #Define sex:CNV
        tmp_df <- df[!grepl(measure,df$X) & grepl("sexM",df$X) & grepl("del",df$X) & df$burden == burden,]
        if (nrow(tmp_df) == 1) {
          out_df[burden,glue("p_del:sexM")] <- tmp_df$P
        }
      } else {
       print(glue("No significant interactions for {file}"))
      }
    }
    summary_df <- rbind(summary_df,out_df)
  }
  summary_df$MATCH <- summary_df$n_sig == summary_df$n_interactions
  summary_df <- summary_df[order(-summary_df$MATCH,summary_df$measure,-summary_df$n_sig),]
  #Write final results
  write.csv(summary_df,file.path("results",pipeline,"stats",pipeline_filter,out_name,"summary_file.csv"))
}
