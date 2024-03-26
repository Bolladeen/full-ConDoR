# %% Load data and plot heatmaps
import yaml
import pandas as pd
import numpy as np
import mosaic.io as mio
import plotly.express as px
from tea.plots import plot_snv_clone
from pathlib import Path

compiled_data_dir = Path("/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_batch2/batch2_data_compiled")

opt_nclones_dir = Path("/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_batch2/condor_downstream/HZ_ete_trees_refined_subclonal_snvs")

output_dir = Path("/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_batch2/condor_downstream/refined_snv_heatmaps")
output_dir.mkdir(parents=True, exist_ok=True)
patient_info_f = "/Users/hzhang/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Iacobuzio_lab/Tapestri_batch2/full-ConDoR/config-MSK.yaml"
with open(patient_info_f, 'r') as f:
    patient_info = yaml.safe_load(f)
patient_names = patient_info["datasets"]

# %%
for patient_name in patient_names:
	if (output_dir / f"{patient_name}-DNA-heatmap.pdf").is_file():
		print(f'[INFO] Skipping {patient_name}.')
		continue

	print(f'[INFO] Processing {patient_name}.')

	pt_h5 = compiled_data_dir / "fillout_h5" / patient_name / f"{patient_name}.patient_wide.genotyped.h5"
	pt = mio.load(pt_h5)

	# Add refined ETE tree clones to H5
	cn_assignment_f = list(opt_nclones_dir.glob(f"{patient_name}/{patient_name}_final*assignment.csv"))[0]
	cn_assignment_df = pd.read_csv(cn_assignment_f, index_col = 0)
	print(f'[INFO] Loaded CN clone assignment file {cn_assignment_f}.')

	# add cn_clone info
	cn_assignment_df['cell_barcode_formatted'] = cn_assignment_df['cell_barcode'] + "-" + cn_assignment_df.index
	cn_assignment_df.set_index('cell_barcode_formatted', inplace=True)
	# for the barcodes in pt that are not in cn_assignment_df, assign them to clone 0
	cn_assignment_df = cn_assignment_df.reindex(pt.dna.barcodes(), fill_value = 0)
	cn_clone_palette = dict(zip(
		np.sort(cn_assignment_df['final_clone_id'].unique()), 
		np.array(px.colors.qualitative.Set3)[np.sort(cn_assignment_df['final_clone_id'].unique())]
	))
	# rename the keys
	cn_clone_palette = {f"CN_clone-{k}": v for k, v in cn_clone_palette.items()}

	pt.dna.row_attrs['label'] = np.array(list(
		map(
			lambda x: f"CN_clone-{int(cn_assignment_df.loc[x, 'final_clone_id'])}", 
			pt.dna.barcodes())
		))
	pt.dna.set_palette(cn_clone_palette)
	pt.cnv.row_attrs['label'] = pt.dna.row_attrs['label']
	pt.cnv.set_palette(cn_clone_palette)
	# pt.cnv.get_gene_names()

	num_cells = pt.dna.shape[0]

	# read in high-quality SNVs
	snv_f = list(compiled_data_dir.glob(f"**/{patient_name}*-voi.hz_curated.txt"))[0]
	snv_df = pd.read_csv(snv_f, sep = '\t', index_col = 0, comment='#')
	snv_ann_map = snv_df['HGVSp'].to_dict()

	# ===== highlight vars with bulk annotation ====='
	# highlight vars
	germline_hom_var_col = '#00cc66' # green
	germline_het_var_col = '#2693ff' # blue
	somatic_var_col = '#ff0000' # red
	likely_artifact_col = '#ffcc00' # yellow

	for var_i in snv_ann_map:
		# germline
		if snv_df.loc[var_i, 'annotation'] == 'germline_HOM':
			snv_ann_map[var_i] = f'<span style="color:{germline_hom_var_col};">' + snv_ann_map[var_i] + '</span>'
		elif snv_df.loc[var_i, 'annotation'] == 'germline_HET':
			snv_ann_map[var_i] = f'<span style="color:{germline_het_var_col};">' + snv_ann_map[var_i] + '</span>'
		elif snv_df.loc[var_i, 'annotation'] == 'bulk_somatic':
			snv_ann_map[var_i] = f'<span style="color:{somatic_var_col};">' + snv_ann_map[var_i] + '</span>'
		elif snv_df.loc[var_i, 'annotation'] == 'likely_artifact':
			snv_ann_map[var_i] = f'<span style="color:{likely_artifact_col};">' + snv_ann_map[var_i] + '</span>'
		else:
			pass

	# plot heatmap
	fig = plot_snv_clone(
		pt,
		sample_name=patient_name,
		story_topic = f'{patient_name}-high_conf_snvs',
		voi = snv_df.index.tolist(),
		attribute = "AF_MISSING",
		ann_map = snv_ann_map
	)
	# # map the x-axis ticks according to snv_ann_map
	# fig.update_xaxes(ticktext = list(map(lambda x: snv_ann_map[x], snv_df.index.tolist())))

	# save the figure
	# write to PDF, with width proportion to the number of SNVs
	fig.write_image(
		str(output_dir / f"{patient_name}-DNA-heatmap.pdf"), 
		width = 500 + 10 * len(snv_df.index.tolist())
	)
	print(f'[INFO] Saved heatmap for {patient_name} to {output_dir / f"{patient_name}-DNA-heatmap.pdf"}.')

# %%
