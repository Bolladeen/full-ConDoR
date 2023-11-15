#!/bin/bash 
#BSUB -J condor_smk
#BSUB -sla IACOBUZC
#BSUB -n 1                # number of core(tasks) for parallel jobs          
#BSUB -R rusage[mem=2]    # expected resorce consumption for memory
#BSUB -W 8:00            # run time limit (hours)
#BSUB -o /home/zhangh5/work/Tapestri_batch2/analysis/condor_pipeline/condor_smk.stdout
#BSUB -eo /home/zhangh5/work/Tapestri_batch2/analysis/condor_pipeline/condor_smk.stderr

if [ -f ~/.bashrc ] ; then
    . ~/.bashrc
fi

mamba activate condor

snakemake -s /home/zhangh5/work/Tapestri_batch2/analysis/full-ConDoR/condor_pipeline.smk \
    --configfile /home/zhangh5/work/Tapestri_batch2/analysis/full-ConDoR/config.yaml \
    --profile lsf \
    --conda-prefix /home/zhangh5/work/Tapestri_batch2/pipeline_results_custom/conda \
    -n
