#!/bin/bash 
#BSUB -J condor_smk
#BSUB -n 1                # number of core(tasks) for parallel jobs          
#BSUB -R rusage[mem=2]    # expected resorce consumption for memory
#BSUB -W 8:00            # run time limit (hours)
#BSUB -o /data/iacobuzc/haochen/Tapestri_batch2/analysis/condor-pipeline/condor_smk.stdout
#BSUB -eo /data/iacobuzc/haochen/Tapestri_batch2/analysis/condor-pipeline/condor_smk.stderr

if [ -f ~/.bashrc ] ; then
    . ~/.bashrc
fi

conda activate condor

snakemake -s /Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_batch2/full-ConDoR/condor_pipeline.smk \
    --configfile /Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_batch2/full-ConDoR/config_MSK_hz_local.yaml \
    --cores 4 \
    --use-conda


