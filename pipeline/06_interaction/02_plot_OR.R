# Libaries
library(tidyverse)
library(glue)
library(stringr)
library(Publish)
library(ggplot2)
library(cowplot)

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
leg <- "none"

orig_names    <- c("Binary_50kb","Binary_100kb","Binary_200kb","N_CNVs","N_genes","Length","LOEUF","pHI","pTS")
publish_names <- c("Presence (50kb)","Presence (100kb)","Presence (200kb)","CNVs","CNV-genes","Length","LOEUF","pHI","pTS")

name_dict <- data.frame(names = publish_names,row.names = orig_names)
print(name_dict)
#For female-only model assign flags that allow for removal of sex terms
sex_flag <- FALSE
sex_t <- " + sex"

input <- "consensus"

analysis_df <- load_data(input,pipeline,pipeline_filter,dx = dx)

if (!(dx == "ksads")) {
  pipeline_filter <- glue("{pipeline_filter}_{dx}")
}

or_cols <- c("Estimate", "Std. Error","CI_0.025","CI_0.975","z value","P","burden")

#Define variable names and measures
burden_types <- orig_names
measures <- c("del","dup")

or_cols <- c("sex","OR","CI_lo","CI_hi","p.value","burden")


for (rm_low_iq in c(TRUE,FALSE)) {
  if (rm_low_iq) {
    out_name <- "burden_IQ"
  } else {
    out_name <- "burden"
  }
  #Create output dirs
  outdir <- file.path("results",pipeline,"stats",pipeline_filter,out_name,disorder,measure)
  dir.create(file.path(outdir,"csvs"), recursive = TRUE,showWarnings = FALSE)

  #Make out summary df
  summary_df <- data.frame(matrix(nrow=length(burden_types)*length(measures)*2,ncol=length(or_cols) + 3)) #nrow =
  colnames(summary_df) <- c("measure","burden_type",or_cols)
  summary_df$measure <- rep(measures,each=length(burden_types)*2)
  summary_df$burden_type <- rep(burden_types,each = 2,length(measures))
  summary_df$sex <- rep(c("male","female"),length(measures)*length(burden_types))

  for (row in 1:nrow(summary_df)) {
    #Load file
    df <- read.csv(glue("{outdir}/csvs/{summary_df[row,'measure']}_sexI_OR.csv"),row.names = 1)
    summary_df[row,or_cols] <- df[df$burden == summary_df[row,"burden_type"] & df$sex == summary_df[row,"sex"],]
  }
  #Remove subset of columns
  summary_df <- summary_df[,!(colnames(summary_df) %in% c("CI.95","burden_type"))]
  #Reassign burden names
  summary_df$burden <- name_dict[summary_df$burden,"names"]

  #Save results
  write.csv(summary_df,glue("{outdir}/summary_OR.csv"))

  plots <- list()
  #Plot results
  n <- 1
  for (result in c("binary","continuous")) {
    print(result)
    if (result == "binary") {
      #Plot binary presence variables
      burden_types_lmt <- c("Presence (50kb)","Presence (100kb)","Presence (200kb)")
      ylim_min <- 0
      ylim_max <- 4
      labs <- c("CNV\n(≥ 50kb)","CNV\n(≥ 100kb)","CNV\n(≥ 200kb)")
      summary_df_tmp <- summary_df[summary_df$burden %in% burden_types_lmt,]
      summary_df_tmp$burden <- factor(summary_df_tmp$burden,levels=burden_types_lmt)
    } else {
      # Plot continuous measures of CNV genes and length
      burden_types_lmt <- c("CNV-genes","Length")
      summary_df_tmp <- summary_df[summary_df$burden %in% burden_types_lmt,]
      summary_df_tmp$burden <- factor(summary_df_tmp$burden,levels=burden_types_lmt)
      ylim_min <- 0.7
      ylim_max <- 1.5
      labs <- c("Number of\nCNV-Genes","Total Length\n(per 100 kb)")
    }
    for (cnv in c("del","dup")) {
      #Select relevant CNV types only
      summary_df_tmp2 <- summary_df_tmp[summary_df_tmp$measure == cnv,]
      if (cnv == "del") {
        lab <- "Deletions"
        leg.pos <- leg
        summary_df_tmp2 <- summary_df_tmp2[summary_df_tmp2$burden != "pTS",]
      } else {
        lab <- "Duplications"
        leg.pos <- leg
        summary_df_tmp2 <- summary_df_tmp2[summary_df_tmp2$burden != "pHI",]
      }
      dodge <- position_dodge(width=0.5)

      #Generate plot
      plots[[n]] <- ggplot(summary_df_tmp2, aes(x = burden, y = OR,color=sex)) +        # ggplot2 plot with confidence intervals
        geom_point(position=dodge) + geom_errorbar(aes(ymin = CI_lo, ymax = CI_hi),position = dodge)  + theme_classic() +
        geom_hline(yintercept=1,linetype=2) + ylim(c(ylim_min,ylim_max)) + labs(x = "",y = "Odds Ratio") + ggtitle(lab) +
        scale_x_discrete(labels = labs) +
        theme(legend.position=leg.pos,axis.text.x = element_text(size = 11),
              axis.title = element_text(size = 11),plot.title = element_text(hjust = 0.5,size = 18))
      n = n + 1
    }
  }
  #Combine plots
  plot_grid(plots[[1]],plots[[2]],plots[[3]],plots[[4]], labels = c("A","","B",""),label_size = 10)
  #Save plots
  ggsave(glue("{outdir}/OR_figure.png"))
}
