#!/bin/bash

module load python
module load R

#Choose filtering criteria
filt=bpd_psych # No bipolar disorder or psychosis

# Run interaction and moderation analyses
for group in consensus 1st_eth; do
  for disorder in ever_ADHD externalizing ODD CD OCD mood anxiety; do
  # for disorder in ever_ADHD mood anxiety conduct OCD; do
    if [[ "$disorder" == "ever_ADHD" ]]; then
      dxs=(ksads cbcl ksad_cbcl)
    else
      dxs=(ksads)
    fi
    #Loop through all diagnostic criteria
    for dx in ${dxs[*]}; do
      #Run logistic regression
      for measure in CNV PRS; do
        echo $ filt $dx $measure $disorder
        Rscript scripts/06_interaction/01_CNV_burden_analysis.R ${group} ${ filt} ${measure} ${disorder} ${dx}
      done
      #Plot ORs
      Rscript scripts/06_interaction/02_plot_OR.R ${group} bpd_psych CNV ${disorder} ${dx}
      #Check for moderation
      if [[ "$disorder" == "ever_ADHD" ]]; then
        Rscript scripts/07_moderation/01_moderation_analysis.R ${group} ${ filt} ${dx}
      fi
    done
  done
done

 filt=bpd_psych
for group in consensus; do
  for disorder in ever_ADHD; do
    for dx in inatt HI; do
      #Run logistic regression
      for measure in CNV; do
        echo $ filt $dx $measure $disorder
        Rscript scripts/06_interaction/01_CNV_burden_analysis.R ${group} ${ filt} ${measure} ${disorder} ${dx}
      done
    done
  done
done


# Run gene enrichment analyses
for group in consensus 1st_eth; do
  for sex in female male all; do
    Rscript scripts/08_bioinformatics/01_gene_enrichment.R consensus ${ filt} ${sex}
  done
done

#Run Mediation Analysis
swarm -f scripts/09_mediation/04_run_mediations.sh -g 120 -t 28 --time=4:00:00 --gres=lscratch:50 --module R --partition=norm,quick
