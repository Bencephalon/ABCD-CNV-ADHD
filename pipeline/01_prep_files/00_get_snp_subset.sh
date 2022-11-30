module load plink


cd /data/NCR_SBRB/jungbt/CNV/ABCD/data/genotype_QCed
#Use the set of SNPs released in the pre-QC'd data as our SNP source for CNV calling
plink --bfile ABCD_release_3.0_QCed --write-snplist
