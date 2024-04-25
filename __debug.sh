# ===== BPA-1 =====
python /data/iacobuzc/haochen/Tapestri_batch2/analysis/full-ConDoR/fast-ConDoR/src/fast-condor.py \
	-i condor_inputs/BPA-1/character_bin_matrix.csv \
	-r condor_inputs/BPA-1/total_readcounts.csv  \
	-v condor_inputs/BPA-1/alt_readcounts.csv \
	-s condor_inputs/BPA-1/germline_mutations.txt \
	-s2 condor_inputs/BPA-1/somatic_mutations.txt  \
	-o condor_outputs/BPA-1/out \
	-m /data/iacobuzc/haochen/Tapestri_batch2/analysis/full-ConDoR/references/amplicon_panel.csv \
	-c /data/iacobuzc/haochen/Tapestri_batch2/batch2_data_BRCA_compiled/SNV_FILLOUT_results \
	-d BPA-1 \
	--scr \
	--cnp /data/iacobuzc/haochen/Tapestri_batch2/batch2_data_BRCA_compiled/falcon_solutions/BPA-1-homdel-nclones=8_solution.amp_clone_profiles.csv \
	--subclonal_mutations /data/iacobuzc/haochen/Tapestri_batch2/analysis/condor-pipeline/manual_subclonal_snvs/BPA-1.subclonal_mutations.yaml 

# # ===== M04 =====
# python /data/iacobuzc/haochen/Tapestri_batch2/analysis/full-ConDoR/fast-ConDoR/src/fast-condor.py \
# 	-i condor_inputs/M04/character_bin_matrix.csv \
# 	-r condor_inputs/M04/total_readcounts.csv  \
# 	-v condor_inputs/M04/alt_readcounts.csv \
# 	-s condor_inputs/M04/germline_mutations.txt \
# 	-s2 condor_inputs/M04/somatic_mutations.txt  \
# 	-o condor_outputs/M04/out \
# 	-m /data/iacobuzc/haochen/Tapestri_batch2/analysis/full-ConDoR/references/amplicon_panel.csv \
# 	-c /data/iacobuzc/haochen/Tapestri_batch2/batch2_data_compiled/SNV_FILLOUT_results \
# 	-d M04 \
# 	--scr \
# 	--cnp /data/iacobuzc/haochen/Tapestri_batch2/analysis/condor-pipeline/0-opt_nclones/M04-homdel-nclones=9_solution.amp_clone_profiles.csv \
# 	--subclonal_mutations /data/iacobuzc/haochen/Tapestri_batch2/analysis/condor-pipeline/manual_subclonal_snvs/M04.subclonal_mutations.yaml 