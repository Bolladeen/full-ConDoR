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

mamba activate condor

snakemake -s /data/iacobuzc/haochen/Tapestri_batch2/analysis/full-ConDoR/condor_pipeline.smk \
    --configfile /data/iacobuzc/haochen/Tapestri_batch2/analysis/full-ConDoR/config-MSK.yaml \
    --use-conda \
    --conda-prefix /data/iacobuzc/haochen/Tapestri_project/TapVarCallSmk/conda \
    --cores all


