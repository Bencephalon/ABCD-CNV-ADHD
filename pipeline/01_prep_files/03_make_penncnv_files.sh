#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=logs/make_penncnv_files.log

#Load Biowulf moduleslogs
module load penncnv
module load R
module load python
module load bedtools

#Set paths
wd=/data/NCR_SBRB/jungbt/CNV/ABCD
genome=/fdb/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa

#Make output dirs
mkdir -p ${wd}/data/model/
cd ${wd}

#Use PennCNV to generate a PFB file from our samples
# compile_pfb.pl ${wd}/data/GSA*/penncnv/USA_Families_*.txt -output ${wd}/data/model/GSA_combined.pfb
echo "Constructing PFB File"

pfb_subs=${wd}/data/model/pfb_subjects.txt
compile_pfb.pl --listfile ${pfb_subs} \
    -output ${wd}/data/model/ABCD.pfb

#Create temporary PFB file with missing header for autodenovo
tail -n +2 ${wd}/data/model/ABCD.pfb > tmp_ABCD.pfb

#Use autodenovo to convert the PFB file into a GCmodel
echo "Constructing GC Model"
bash ~/autodenovo/penncnv_create_gcmodel.sh \
                          ${wd}/tmp_ABCD.pfb \
                          ${genome} \
                          ${wd}/data/model/ABCD.hg19.gcmodel
#Remove temp PFB
rm tmp_ABCD.pfb
