#!/usr/bin/env bash


# April 2022
# Kwangmi An
# for WNH
# Modified BJ - May 2022

module load R
module load plink/2.3-alpha

wd=/data/NCR_SBRB/jungbt/CNV/ABCD
keep_file=${wd}/data/PLINK/keep_subjects_pca.WNH.txt
out_prs=${wd}/results/PRS/
out_pca=${wd}/results/PCA/WNH/
plink_file=${wd}/data/PLINK/ABCD_release_3.0_QCed.for_Ben

mkdir -p ${out_prs}
mkdir -p ${out_pca}

#Run PRS analysis
Rscript /data/NCR_SBRB/software/PRSice_2.3.3/PRSice.R \
--prsice /data/NCR_SBRB/software/PRSice_2.3.3/PRSice_linux \
--base /data/NCR_SBRB/Ahn/ABCD/PRS/adhd.2019 \
--base-info INFO:0.9 \
--target ${plink_file} \
--keep ${keep_file} \
--all-score \
--binary-target F \
--stat OR \
--print-snp T \
--a1 A1 \
--a2 A2 \
--bp BP \
--chr CHR \
--pvalue P \
--clump-p 1 \
--clump-r2 0.1 \
--clump-kb 250 \
--fastscore T \
--bar-levels 1,0.5,0.1,0.2,0.05,0.01,0.001,0.0001,0.00001,0.00001,0.00000005 \
--no-regress \
--thread 28 \
--out ${out_prs}/ABCD_release_3.0_QCedfor_Ben.PRS_adhd2019

#Run PCA Calculation - WNH
plink2 --bfile ${out_pca}/../whole_sample/snp_for_ancestry \
       --keep ${keep_file} \
       --maf 0.10 --autosome  --pca \
       --out ${out_pca}/pca_analysis.WNH
