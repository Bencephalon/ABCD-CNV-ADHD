# swarm -f scripts/run_mediations.sh -g 120 -t 28 --time=4:00:00 --gres=lscratch:50 --module R --partition=norm,quick

#Test Cognitive Mediators
Rscript scripts/09_mediation/02_cog_mediation.R consensus bpd_psych female
Rscript scripts/09_mediation/02_cog_mediation.R consensus bpd_psych
Rscript scripts/09_mediation/02_cog_mediation.R consensus bpd_psych male
Rscript scripts/09_mediation/02_cog_mediation.R 1st_eth   bpd_psych female
Rscript scripts/09_mediation/02_cog_mediation.R 1st_eth   bpd_psych
Rscript scripts/09_mediation/02_cog_mediation.R 1st_eth   bpd_psych male
#Test Neural Mediators
Rscript scripts/09_mediation/03_mri_mediation.R consensus bpd_psych female targeted
Rscript scripts/09_mediation/03_mri_mediation.R consensus bpd_psych male   targeted
Rscript scripts/09_mediation/03_mri_mediation.R consensus bpd_psych all    targeted
Rscript scripts/09_mediation/03_mri_mediation.R 1st_eth   bpd_psych female targeted
Rscript scripts/09_mediation/03_mri_mediation.R 1st_eth   bpd_psych male   targeted
Rscript scripts/09_mediation/03_mri_mediation.R 1st_eth   bpd_psych all    targeted
