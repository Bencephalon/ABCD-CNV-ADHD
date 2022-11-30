#!/bin/bash
#SBATCH -n 24
#SBATCH --mem=120G
#SBATCH --time=4:00:00
#SBATCH --partition=norm,quick
#SBATCH --output=logs/penncnv_postprocess.log
#Load Biowulf modules
module load penncnv

#Specify paths
wd=/data/NCR_SBRB/jungbt/CNV/ABCD
indir=${wd}/results/penncnv/
outdir=${wd}/results/penncnv/merged/
pfb=${wd}/data/model/ABCD.pfb

#Set QC thresholds
qcbafdrift=0.01
qcwf=0.05
qclrrsd=0.35
frac=0.5

#Make outdir
rm -r ${outdir}
mkdir -p ${outdir}

#Merge all swarm commands into a single CNV file
cat ${indir}/*_ABCD.log >> ${outdir}/ABCD.log
cat ${indir}/*_ABCD.rawcnv >> ${outdir}/ABCD.rawcnv

#Use PennCNV to filter out low quality CNVs
filter_cnv.pl ${outdir}/ABCD.rawcnv \
      --qcbafdrift ${qcbafdrift} -qclogfile ${outdir}/ABCD.log \
      -qcwf ${qcwf} -qclrrsd ${qclrrsd} \
      -out ${outdir}/ABCD.filtered.rawcnv \
      -qcsumout ${outdir}/ABCD.filtered.qcsum
