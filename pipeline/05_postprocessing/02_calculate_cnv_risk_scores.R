
# if (!"BiocManager" %in% rownames(installed.packages()))
#      install.packages("BiocManager")
# BiocManager::install(c("AnnotationHub", "Homo.sapiens","TxDb"))
library(glue)
library(AnnotationHub)
library(GenomicRanges)
library(IRanges)
library(Repitools)
library(stringr)

#Load annotation hub
ah <- AnnotationHub()
print(ah[ah$rdataclass == "GRanges" & ah$dataprovider == "GENCODE" & ah$genome == "GRCh37",])
grs <- ah["AH75189",]
genes <- grs[[1]]
#Remove entries with no gene name
genes <- genes[!is.na(genes$Gene),]

#Make gene list for plink
genes_df <- annoGR2DF(genes)
#Chr (without "chr"), start, stop, name
out_genes_df <- data.frame("CHR" = unlist(lapply(genes_df[,1],str_replace,"chr","")),"START" = genes_df[,2],"END" = genes_df[,3],"NAME" = genes_df[,10])
write.table(out_genes_df,glue("/data/NCR_SBRB/jungbt/CNV/hg19_regions/gene_lists/PLINK_gene_table.genes"),sep = "\t",quote = FALSE,row.names = FALSE, col.names = FALSE)

#Load CNVs
df <- read.csv("results/consensus/processed_cnvs_rare.csv")
df$Chr <- paste("chr",df$Chr,sep="")
df$Length <- df$End - df$Start
df$coordcnv <- glue("{df$Chr}:{df$Start}-{df$End}")

#Create a list of genes overlapping CNVs
cnv = GenomicRanges::makeGRangesFromDataFrame(df,keep.extra.columns=TRUE)
olaps = GenomicRanges::findOverlaps(cnv, genes)
long_annotated = cnv[queryHits(olaps)]
long_annotated$gene_id = genes[subjectHits(olaps)]$Gene

#Load pHI and pTS scores
dosage_sensitivity_scores <- read.table("data/other/Collins_rCNV_2022.dosage_sensitivity_scores.tsv")
colnames(dosage_sensitivity_scores) <- c("Gene","pHaplo","pTriplo")
rownames(dosage_sensitivity_scores) <- dosage_sensitivity_scores$Gene
#Assign pHI/pTS scores
long_annotated$pHI   <- dosage_sensitivity_scores[as.matrix(long_annotated$gene_id),"pHaplo"]
long_annotated$pTS   <- dosage_sensitivity_scores[as.matrix(long_annotated$gene_id),"pTriplo"]


#Load LOEUF scores
loeuf_df <- read.table("/data//NCR_SBRB/jungbt/CNV/hg19_regions/gnomad.v2.1.1.lof_metrics.by_gene.txt",sep="\t",header=TRUE)
loeuf_df <- loeuf_df[,c("gene","oe_lof_upper")]
loeuf_df <- aggregate(loeuf_df$oe_lof_upper,by = list(loeuf_df$gene), data = loeuf_df, FUN = mean)
colnames(loeuf_df) <- c("gene","LOEUF")
rownames(loeuf_df) <- loeuf_df$gene
#Assign LOEUF scores
long_annotated$LOEUF <- loeuf_df[as.matrix(long_annotated$gene_id),"LOEUF"]

#Save gene-based info to new var
cnv_genes <- as.data.frame(long_annotated)
#Calculate the inverse of the LOEUF scores
cnv_genes$iLOEUF <- 1/cnv_genes$LOEUF
#Load in subject df
sub_df <- read.csv("results/consensus/subject_info.csv")

#Calculate number of genes
out_sample <- data.frame(matrix(nrow=nrow(sub_df),ncol = 0))
rownames(out_sample) <- sub_df$ID

dir.create("results/consensus/analysis_tables/genes/",recursive = TRUE)

#For each CNV type
for (cnv_type in c("del","dup","both")) {
  print(cnv_type)
  #Subset CNV data
  if (cnv_type == "del") {
    cnvs_tmp <- df[df$CNV_Value < 2,]
    cnv_genes_tmp <- cnv_genes[cnv_genes$CNV_Value < 2,]
  } else if (cnv_type == "dup") {
    cnvs_tmp <- df[df$CNV_Value > 2,]
    cnv_genes_tmp <- cnv_genes[cnv_genes$CNV_Value > 2,]
  } else {
    cnvs_tmp <- df
    cnv_genes_tmp <- cnv_genes
  }
  # Calculate measures of CNV burden/presence
  #Binary_50kb
  cnvs_tmp_thresh <- cnvs_tmp[cnvs_tmp$Length >= 50000,]
  out_sample[,glue("Binary_50kb_{cnv_type}")] <- rownames(out_sample) %in% cnvs_tmp_thresh$Sample_ID
  #Binary_100kb
  cnvs_tmp_thresh <- cnvs_tmp[cnvs_tmp$Length >= 100000,]
  out_sample[,glue("Binary_100kb_{cnv_type}")] <- rownames(out_sample) %in% cnvs_tmp_thresh$Sample_ID
  #Binary_200kb
  cnvs_tmp_thresh <- cnvs_tmp[cnvs_tmp$Length >= 200000,]
  out_sample[,glue("Binary_200kb_{cnv_type}")] <- rownames(out_sample) %in% cnvs_tmp_thresh$Sample_ID
  for (sample in rownames(out_sample)) {
    # print(n)
    #1) Calculate based on CNVs (not CNV-genes)
    #N_CNVs
    out_sample[sample,glue("N_CNVs_{cnv_type}")]  <- nrow(cnvs_tmp[cnvs_tmp$Sample_ID == sample,])
    #Length
    out_sample[sample,glue("Length_{cnv_type}")]  <- sum(cnvs_tmp[cnvs_tmp$Sample_ID == sample,"Length"],na.rm=TRUE)/100000
    #2) Calculate based on CNV-genes
    #N_genes
    out_sample[sample,glue("N_genes_{cnv_type}")] <- nrow(cnv_genes_tmp[cnv_genes_tmp$Sample_ID == sample,])
    #pHI
    out_sample[sample,glue("pHI_{cnv_type}")]     <- sum(cnv_genes_tmp[cnv_genes_tmp$Sample_ID == sample,"pHI"],na.rm=TRUE)
    #pTS
    out_sample[sample,glue("pTS_{cnv_type}")]     <- sum(cnv_genes_tmp[cnv_genes_tmp$Sample_ID == sample,"pTS"],na.rm=TRUE)
    #LOEUF
    out_sample[sample,glue("LOEUF_{cnv_type}")]     <- sum(cnv_genes_tmp[cnv_genes_tmp$Sample_ID == sample,"iLOEUF"],na.rm=TRUE)
  }
  #Make gene table
  cnv_freq_table <- as.data.frame.matrix(table(cnv_genes_tmp$Sample_ID,cnv_genes_tmp$gene_id))
  samples <- rownames(cnv_freq_table)

  # Collapse identical CNV-genes
  cnv_freq_table <- setNames(as.data.frame(unique(as.list(cnv_freq_table))),
           sapply(split(names(cnv_freq_table), match(as.list(cnv_freq_table), unique(as.list(cnv_freq_table)))),
                  paste, collapse="+"))
  rownames(cnv_freq_table) <- samples
  #Add in subjects without CNVs
  cnv_freq_table[rownames(out_sample)[!(rownames(out_sample) %in% samples)],] <- 0
  #Binarize
  cnv_freq_table <- cnv_freq_table > 0

  #Save gene table
  write.csv(cnv_freq_table,glue("results/consensus/analysis_tables/genes/gene_table_{cnv_type}.csv"))
}
dir.create("results/consensus/analysis_tables/burden/",recursive = TRUE)
#Write burden table
write.csv(out_sample,"results/consensus/analysis_tables/burden/burden_table.csv")
