# %% Load data and plot heatmaps
import yaml
import pandas as pd
import numpy as np
import mosaic.io as mio
import plotly.express as px
from tea.plots import plot_snv_clone
from tea.utils import sort_for_var
from pathlib import Path

compiled_data_dir = Path("/lila/data/iacobuzc/haochen/Tapestri_batch2/batch2_data_BRCA_compiled")
opt_nclones_dir = Path("/lila/data/iacobuzc/haochen/Tapestri_batch2/analysis/condor_downstream/HZ_ete_trees_refined_subclonal_snvs")

output_dir = Path("/lila/data/iacobuzc/haochen/Tapestri_batch2/analysis/condor_downstream/refined_snv_heatmaps")
output_dir.mkdir(parents=True, exist_ok=True)
patient_info_f = "/lila/data/iacobuzc/haochen/Tapestri_batch2/analysis/full-ConDoR/config-MSK.yaml"
with open(patient_info_f, 'r') as f:
    patient_info = yaml.safe_load(f)
patient_names = patient_info["datasets"]

# @HZ 2024-03-27
# highlight the gBRCAs
VARS_TO_HIGHLIGHT = ["chr13:32936733:A/T", "chr13:32914437:GT/G", "chr13:32914288:TAACC/T"]

# get HOMDEL
amplicon_params = Path('/lila/data/iacobuzc/haochen/Tapestri_project/tap_cn_calling/train-normals/train-combined_8_normals/NB_train-combined_8_normals-results.gene_added.csv')
amp_params_df = pd.read_csv(
    amplicon_params,
    index_col= 0,
)
# %%
for patient_name in patient_names:
	if (output_dir / f"{patient_name}-DNA-heatmap.pdf").is_file():
		print(f'[INFO] Skipping {patient_name}.')
		continue

	print(f'[INFO] Processing {patient_name}.')

	pt_h5 = (compiled_data_dir / "SNV_FILLOUT_results" / patient_name).glob(f"{patient_name}.*.h5")
	pt_h5 = list(pt_h5)[0]
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
	BARCODE="barcode"
	pt.cnv.add_row_attr(BARCODE, pt.dna.barcodes())
	pt.cnv.row_attrs['label'] = pt.dna.row_attrs['label']
	pt.cnv.set_palette(cn_clone_palette)
	pt.cnv.get_gene_names(amplicon_params, gene_name_col = 'gene')
	pt.cnv.var = pd.DataFrame.from_dict(pt.cnv.col_attrs).set_index("id")
 
	raw_rc = pt.cnv.get_attribute('read_counts', constraint='row')
	normalized_rc = raw_rc / amp_params_df.loc[raw_rc.columns]['amplicon_factor'] / raw_rc.sum(axis=1).values[:, None]
	# add to sample
	pt.cnv.add_layer('normalized_read_counts', normalized_rc.values)
	normalized_rc_binary = (normalized_rc == 0).astype(int)
	pt.cnv.add_layer('normalized_read_counts_binary[zero/nonzero]', normalized_rc_binary.values)

	num_cells = pt.dna.shape[0]

	# read in high-quality SNVs
	snv_f = list(compiled_data_dir.glob(f"**/{patient_name}*-voi.hz_curated.txt"))[0]
	snv_df = pd.read_csv(snv_f, sep = '\t', index_col = 0, comment='#')
	# get rid of the artifact vars, which has "artifact" in their annotation column
	snv_df = snv_df[~snv_df["annotation"].str.contains("artifact", na=False)]
	# then, the vars with no annotation would be "novel"
	snv_df.loc[snv_df["annotation"].isna(), "annotation"] = "novel"
	snv_ann_map = snv_df['HGVSp'].to_dict()

	# ===== highlight vars with bulk annotation ====='
	# highlight vars
	germline_hom_var_col = '#00cc66' # green
	germline_het_var_col = '#2693ff' # blue
	somatic_var_col = '#ff0000' # red
	# likely_artifact_col = '#ffcc00' # yellow
	highlight_color = '#ffcc00'

	for var_i in snv_ann_map:
		if var_i in VARS_TO_HIGHLIGHT:
			snv_ann_map[var_i] = f'<span style="color:{highlight_color};">' + snv_ann_map[var_i] + '</span>'
		# germline
		elif snv_df.loc[var_i, 'annotation'] == 'germline_HOM':
			snv_ann_map[var_i] = f'<span style="color:{germline_hom_var_col};">' + snv_ann_map[var_i] + '</span>'
		elif snv_df.loc[var_i, 'annotation'] == 'germline_HET':
			snv_ann_map[var_i] = f'<span style="color:{germline_het_var_col};">' + snv_ann_map[var_i] + '</span>'
		elif "somatic" in snv_df.loc[var_i, 'annotation']:
			snv_ann_map[var_i] = f'<span style="color:{somatic_var_col};">' + snv_ann_map[var_i] + '</span>'
		# elif "artifact" in snv_df.loc[var_i, 'annotation']:
		# 	snv_ann_map[var_i] = f'<span style="color:{likely_artifact_col};">' + snv_ann_map[var_i] + '</span>'
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
	fig.write_image(
		str(output_dir / f"{patient_name}-DNA-heatmap.png"), 
		width = 500 + 10 * len(snv_df.index.tolist()),
  		scale = 3,
	)
	print(f'[INFO] Saved heatmap for {patient_name} to {output_dir / f"{patient_name}-DNA-heatmap.pdf"}.')

	attribute="AF_MISSING"
	# also make sure to retrieve the sorted_bars 
	del pt.dna.__dict__["_Assay__heatmap"]
	sorted_bars = sort_for_var(
		dna = pt.dna,
		vars = snv_df.index.tolist(),
		attribute = attribute,
		method = "hier"
		)

	sc_cnv_heatmap = pt.cnv.heatmap(
		'normalized_read_counts_binary[zero/nonzero]',
		features = ['CDKN2A', 'KRAS', 'BRCA2', 'TP53', 'SMAD4'],
		bars_order = sorted_bars,
	)

	# update colorscale title to be: "homdel"
	# update colorscale to be two discrete values: 0, 1
	# decrease the height of the colorscale
	sc_cnv_heatmap.layout.coloraxis.colorbar.title = 'zero-reads'
	sc_cnv_heatmap.layout.coloraxis.colorbar.len = 0.2
	sc_cnv_heatmap.layout.coloraxis.colorbar.tickvals = [0, 1]


	sc_cnv_heatmap.layout.coloraxis.colorscale = [
		[0, 'rgb(0, 0, 0)'], 
		[1, 'rgb(255, 225, 0)']
	]
	sc_cnv_heatmap.update_layout(font_family="Arial")
	sc_cnv_heatmap.show()
	sc_cnv_heatmap.write_image(
		output_dir / f'{patient_name}-cnv_heatmap-binarized_select_genes.png',
		height=800, width=800, scale=3,
		)


# %%
