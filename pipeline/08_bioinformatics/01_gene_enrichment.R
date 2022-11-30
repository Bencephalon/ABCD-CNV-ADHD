library(tidyverse)
library(readxl)
library(glue)

source("scripts/99_util/util_load_cnvs.R")

set.seed(2514)

#Load pipeline
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  pipeline = "consensus"
  pipeline_filter = "bpd_psych"
  sex = "all"
  dx = "ksads"
} else if (length(args) == 1) {
  pipeline = args[1]
  pipeline_filter = "bpd_psych"
  sex = "all"
  dx = "ksads"
} else if (length(args) == 2) {
  pipeline = args[1]
  pipeline_filter = args[2]
  sex = "all"
  dx = "ksads"
} else if (length(args) == 3) {
  pipeline = args[1]
  pipeline_filter = args[2]
  sex = args[3]
  dx = "ksads"
} else if (length(args) == 4) {
  pipeline = args[1]
  pipeline_filter = args[2]
  sex = args[3]
  dx = args[4]
}

print(pipeline)
print(pipeline_filter)

input <- "consensus"

#Load in analysis file
analysis_df <- load_data(input,pipeline,pipeline_filter,sex = sex,dx = dx)

if (!(sex == "all")) {
  pipeline_filter <- glue("{pipeline_filter}_{sex}")
}

run_model <- function(df,loci,event,out,burden_any = FALSE) {
  #Store p_values
  r_names = loci
  c_names = c("p","OR","CI_low","CI_high")
  p <- data.frame(matrix(0,ncol = length(c_names), nrow = length(r_names)))
  rownames(p) <- r_names
  colnames(p) <- c_names
  c_names = c("case","control","total","prop_case","prop_control","prop_total")
  freq_df <- data.frame(matrix(0,ncol = length(c_names), nrow = length(r_names)))
  rownames(freq_df) <- r_names
  colnames(freq_df) <- c_names
  # print(freq_df)
  if (length(r_names) == 0) {
    print("No CNVs overlap with given loci. Skipping")
    return()
  }
  # Loop through each locus
  if (burden_any & length(loci) > 1) {
    #Calculate compbined metrics
    df$SUBSET <- rowSums(df[,loci], na.rm = TRUE) > 0
    df$ALL <- df[,glue("N_genes_{event}")] > 0
    df$EXCLUDE <- df$ALL & (df$SUBSET == FALSE)
    loci <- c(loci, "SUBSET","EXCLUDE","ALL")
  }
  for (locus in loci) {
    # Create Carrier Variable for each subject
    df$carrier_status <- df[[locus]]
    # Generate descriptive stats
    freq_df[locus,"case"]    <- sum(df$carrier_status & (df$ever_ADHD == "ADHD"))
    freq_df[locus,"control"] <- sum(df$carrier_status & (df$ever_ADHD == "non-ADHD"))
    freq_df[locus,"total"]   <- sum(df$carrier_status)
    freq_df[locus,"prop_case"]    <- freq_df[locus,"case"]/sum(df$ever_ADHD == "ADHD")
    freq_df[locus,"prop_control"] <- freq_df[locus,"control"]/sum(df$ever_ADHD == "non-ADHD")
    freq_df[locus,"prop_total"]   <- freq_df[locus,"total"]/dim(df)[1]
    # Creat contingency table
    cont.table <- table(df$carrier_status, df$ever_ADHD)
    if (nrow(cont.table) > 1) {
      #Run fisher's exact test
      model <- fisher.test(cont.table)
      # Store P-Values, OR and CI
      p[locus,"p"] <- model$p.value
      p[locus,"OR"] <-  model$estimate
      p[locus,"CI_low"] <-  model$conf.int[1]
      p[locus,"CI_high"] <-  model$conf.int[2]
    } else {
      #Don't run fisher's exact test
      p[locus,c("p","OR","CI_low","CI_high")] <- NA
    }
  }
  # Multiple Correction
  freq_df$OR <- p[,"OR"]
  freq_df$p <- p[,"p"]
  freq_df$CI_low <- p[,"CI_low"]
  freq_df$CI_high <- p[,"CI_high"]
  freq_df$p_FDR <- p.adjust(freq_df$p,"BH")
  freq_df$p_FWER <- p.adjust(freq_df$p,"holm")
  freq_df <- freq_df[order(freq_df$p),]

  #Write model results
  write.csv(freq_df,glue("{out}.csv"))
  #Print summary
  print(glue("Loci: {dim(freq_df)[1]}"))
  print(glue("Marginally Significant Loci (p < 0.1): {nrow(freq_df[freq_df$p <= 0.1,])}"))
  print(freq_df[freq_df$p <= 0.1,c("prop_case","prop_control","prop_total","OR","p")])

  return(freq_df)
}

define_gene_list <- function(df,file) {
  #Load existing gene list
  gene_list <- read.csv(file,header=FALSE,col.names=c("gene"))
  gene_list$match = NA
  #Loop through each gene in list
  for (n in 1:nrow(gene_list)) {
    gene <- gene_list$gene[n]
    #Find matching gene in gene table
    x <- grepl(glue("\\b{gene}\\b"), colnames(gene_df))
    #Add actual gene name to gene list
    gene_list[n,"match"] <- colnames(gene_df)[x]
  }
  #Remove genes with no matches
  gene_list <- gene_list[!is.na(gene_list$match),"match"]
  return(unique(gene_list))
}

#Specify output directory name
outdir <- file.path("results",pipeline,"stats",pipeline_filter,"gwas")
dir.create(glue("{outdir}/csvs/"), recursive = TRUE,showWarnings = FALSE)

for (event in c("del","dup","both")) {
  print(event)
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  #Attach genes
  gene_df <- read.csv(glue("results/{input}/analysis_tables/genes/gene_table_{event}.csv"),row.names = 1,check.names = FALSE)
  analysis_df_tmp <- analysis_df
  #Subset gene_df
  gene_df <- gene_df[rownames(gene_df) %in% rownames(analysis_df_tmp),]
  gene_df <- gene_df[,colSums(gene_df) > 0]
  analysis_df_tmp[,colnames(gene_df)] <- gene_df[rownames(analysis_df_tmp),colnames(gene_df)]

  # Filter out very low frequency CNVs
  gene_freq <- colSums(gene_df)
  t <- 1000
  gene_freq <- gene_freq[gene_freq >= nrow(analysis_df_tmp)/t]

  out_pref <- glue("{outdir}/csvs/{event}_gwas")
  freq_df <- run_model(analysis_df_tmp,names(gene_freq),event,out_pref,burden_any = FALSE)

  print("===================================================================")

  gene_list_combined <- c()
  #Load gene lists
  gene_list_harich        <- define_gene_list(gene_df,"/data/NCR_SBRB/jungbt/CNV/hg19_regions/gene_lists/gene_list_harich.txt")
  gene_list_demontis_2022_finemap <- define_gene_list(gene_df,"/data/NCR_SBRB/jungbt/CNV/hg19_regions/gene_lists/gene_list_demontis_2022_finemap.txt")
  gene_list_combined      <- unique(c(gene_list_harich,gene_list_demontis_2022_finemap))

  gene_lists <- list(gene_list_harich,gene_list_demontis_2022_finemap)
  out_names <- c(glue("{event}_harich"),glue("{event}_demontis_2022_finemap"),glue("{event}_combined"),glue("{event}_GO0098793"))

  #Reload gene df
  gene_df <- read.csv(glue("results/{input}/analysis_tables/genes/gene_table_{event}.csv"),row.names = 1,check.names = FALSE)
  gene_df <- gene_df[rownames(gene_df) %in% rownames(analysis_df_tmp),]

  for(n in 1:length(gene_lists)) {
    gene_list <- gene_lists[[n]]
    if (length(gene_list) > 0) {
      print("    running...")
      #Create analysis dataframe
      analysis_df_tmp <- analysis_df
      analysis_df_tmp[,gene_list] <- gene_df[rownames(analysis_df_tmp),gene_list]
      #Specify output name
      out_pref <- glue("{outdir}/csvs/{out_names[n]}")
      #Run exact tests
      freq_df <- run_model(analysis_df_tmp,gene_list,event,out_pref,burden_any = TRUE)
    }
  }
}
