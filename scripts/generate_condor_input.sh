mamba activate mosaic-custom

script=/home/zhangh5/work/Tapestri_batch2/analysis/condor_pipeline/scripts/generate_condor_input.py

python ${script} \
    -d M04 \
    -l /home/zhangh5/work/Tapestri_batch2/batch2_data_compiled/SNV_FILLOUT_results/ \
    -i /home/zhangh5/work/Tapestri_batch2/analysis/condor_pipeline/opt_nclones/M04-homdel-nclones=9_solution.cell_assignments.csv \
    -snvs /home/zhangh5/work/Tapestri_batch2/analysis/snv-cnv-combined/M04-patient/sc_heatmaps/all_vars/mut_prev=0.01/all_vars-voi.hz.txt \
    -v /home/zhangh5/work/Tapestri_batch2/analysis/condor_pipeline/M04/condor_inputs/character_vaf_matrix.csv \
    -m /home/zhangh5/work/Tapestri_batch2/analysis/condor_pipeline/M04/condor_inputs/character_bin_matrix.csv \
    -a /home/zhangh5/work/Tapestri_batch2/analysis/condor_pipeline/M04/condor_inputs/alt_readcounts.csv \
    -t /home/zhangh5/work/Tapestri_batch2/analysis/condor_pipeline/M04/condor_inputs/total_readcounts.csv \
    -g /home/zhangh5/work/Tapestri_batch2/analysis/condor_pipeline/M04/condor_inputs/germline_mutations.txt \
    -s /home/zhangh5/work/Tapestri_batch2/analysis/condor_pipeline/M04/condor_inputs/somatic_mutations.txt