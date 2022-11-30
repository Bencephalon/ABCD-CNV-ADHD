#!/bin/bash

module load python
module load R

# NOTE: Need script for SNP level filtering

####################### Step 1: Format and Prepare Data #######################
# Extract only the SNPs that survived QC metrics for CNV detection
bash    scripts/01_prep_files/00_get_snp_subset.sh
# Reformat ABCD probe data into PennCNV-compatible format
python  scripts/01_prep_files/01_reformat_abcd_cnv_data.py
# Collect demographic and phenotypic information for each subject. Define ADHD
# and other NDDs. Generate PFB file for PennCNV
python  scripts/01_prep_files/02_make_pheno_file.py
# Generate PFB and GCmodel files for PennCNV
bash    scripts/01_prep_files/03_make_penncnv_files.sh
# Generate swarm commands to run all data through PennCNV
python  scripts/01_prep_files/04_make_penncnv_swarm.py
# Reformat data for QuantiSNP. Generate swarm commands to run all data
# through QuantiSNP
python  scripts/01_prep_files/05_make_QT_calls.py

########################## Step 2: Run CNV Detection ##########################
#Run CNV calling through PennCNV
bash    scripts/02_CNV_detection/01_run_penncnv.sh
#Run CNV calling through QuantiSNP
bash    scripts/02_CNV_detection/02_run_QT.sh

########################## Step 3: Process CNV Calls ##########################
# Filter PennCNV outputs for quality (part 1)
bash    scripts/03_CNV_processing/01_postprocess_penncnv.sh
# Filter PennCNV outputs for quality (part 2)
python  scripts/03_CNV_processing/02_postprocess_penncnv.py
# Filter QuantiSNP outputs for quality and reformat outputs
python  scripts/03_CNV_processing/03_process_QT.py
#Filter our list of CNVs by size, consensus, location, and rarity
python  scripts/03_CNV_processing/04_filter_cnvs.py consensus
#Further filter CNVs by size
Rscript scripts/03_CNV_processing/05_filter_rare_cnvs.R consensus

########################## Step 4: Calculate PRS ##########################
# Select one member per family for PCA analysis
Rscript scripts/04_PRS_detection/01_abcd_keep_list.R
# Calculate PCA and run ADXMIXTURE calculations
bash    scripts/04_PRS_detection/02_snp_pcs.sh
# Generate a list of white, non-hispanic individuals for PRS calculations
Rscript scripts/04_PRS_detection/03_abcd_keep_list_WNH.R
# Calculate PCA and PRS for white,non-Hispanic subcohort
bash scripts/04_PRS_detection/04_PCA_PRS_wnh.sh

############# Step 5: Incorporate and Prepare Data Prior to Analysis #############
# Incorporate the ancestry PCs into the subject data
python  scripts/05_postprocessing/01_incorporate_pcs.py
#Calculate measures of CNV burden and presence for each sample
Rscript scripts/05_postprocessing/02_calculate_cnv_risk_scores.R

#Prepare MRI features for mediation analysis
Rscript scripts/09_mediation/01_prepare_mri.R

########################## Step 3: Analysis ##########################
Rscript scripts/run_analyses.sh
