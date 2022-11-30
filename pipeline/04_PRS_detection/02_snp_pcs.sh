#!/bin/bash
#SBATCH -n 32
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/make_snp_pc_files.txt

module load plink/2.3-alpha
module load admixture
wd=/data/NCR_SBRB/jungbt/CNV/ABCD
out_pca=${wd}/results/PCA/whole_sample/

mkdir -p ${out_pca}

plink_file=/data/NCR_SBRB/jungbt/CNV/ABCD/data/PLINK/ABCD_release_3.0_QCed.for_Ben
keep_file=/data/NCR_SBRB/jungbt/CNV/ABCD/data/PLINK/keep_subjects_pca.txt


module load plink

#Edit data set to make all singletons founders (we have previously selected one-per-family)
plink --bfile ${plink_file} \
       --keep ${keep_file} --make-founders \
       --make-bed --out ${out_pca}/founders

module load plink/2.3-alpha

#LD pruning
plink2 --bfile ${out_pca}/founders --keep  ${keep_file} \
       --maf 0.10 --autosome --indep-pairwise 50 0.2 \
       --out ${out_pca}/pruned_snps

#Write out final list of SNPs and subjects
plink2 --bfile ${out_pca}/founders --keep  ${keep_file} \
       --maf 0.10 --autosome \
       --exclude ${out_pca}/pruned_snps.prune.out \
       --make-bed --out ${out_pca}/snp_for_ancestry --write-snplist

#PCA Calculation
plink2 --bfile ${out_pca}/snp_for_ancestry \
       --keep ${keep_file} \
       --maf 0.10 --autosome  --pca \
       --out ${out_pca}/pca_analysis

cd ${out_pca}

rm admixture.swarm

#Run Admixture calculations using 2 - 8 clusters
for K in 6 8 7 5 4 3 2; do
  #ADMIXTURE Population calculations
  echo "admixture -j36 -C=0.1 --cv snp_for_ancestry.bed $K | tee log${K}.out" >> admixture.swarm
done

swarm -f admixture.swarm -m admixture -g 32 -t 36 --time=36:00:00

# Print out cross-validation results
# grep -h CV log*.out
