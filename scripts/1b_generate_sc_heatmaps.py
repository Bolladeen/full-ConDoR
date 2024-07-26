# %% IMPORT MODULES
from pathlib import Path
import pandas as pd
import numpy as np
import mosaic.io as mio
import plotly.express as px
from tea.plots import plot_snv_clone
import argparse
import sys

def main(args):
    patient_name = args.patient_name
    pt_h5 = args.patient_h5
    snv_f = args.snv_f
    cn_assignment_f = args.clone_assignment
    output_dir = Path(args.output_dir)

    pt = mio.load(pt_h5)

    # Add raw FALCON results to H5
    cn_assignment_df = pd.read_csv(cn_assignment_f, index_col = 0)
    print(f'[INFO] Loaded CN clone assignment file {cn_assignment_f}.')

    if args.post_condor:
        cn_assignment_df['clone_id'] = cn_assignment_df["final_clone_id"]

    # add cn_clone info
    cn_assignment_df['cell_barcode_formatted'] = cn_assignment_df['cell_barcode'] + "-" + cn_assignment_df.index
    cn_assignment_df.set_index('cell_barcode_formatted', inplace=True)

    if args.post_condor:
        # also fill in the diploid cells
        # assign cells from pt that are not in cn_assignment_df to clone 0
        cn_assignment_df = cn_assignment_df.reindex(pt.dna.barcodes(), fill_value=0)

    cn_clone_palette = dict(zip(
        np.sort(cn_assignment_df['clone_id'].unique()), 
        np.array(px.colors.qualitative.Set3)[np.sort(cn_assignment_df['clone_id'].unique())]
        ))
    # rename the keys
    cn_clone_palette = {f"CN_clone-{k}": v for k, v in cn_clone_palette.items()}

    pt.dna.row_attrs['label'] = np.array(list(
        map(
            lambda x: f"CN_clone-{int(cn_assignment_df.loc[x, 'clone_id'])}", 
            pt.dna.barcodes())
        ))
    pt.dna.set_palette(cn_clone_palette)
    pt.cnv.row_attrs['label'] = pt.dna.row_attrs['label']
    pt.cnv.set_palette(cn_clone_palette)
    # pt.cnv.get_gene_names()


    # read in high-quality SNVs
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
        story_topic = 'high_conf_snvs_pre_condor_clone_assignment',
        voi = snv_df.index.tolist(),
        attribute = "AF_MISSING",
        ann_map = snv_ann_map
    )
    # # map the x-axis ticks according to snv_ann_map
    # fig.update_xaxes(ticktext = list(map(lambda x: snv_ann_map[x], snv_df.index.tolist())))

    # save the figure
    # write to PDF, with width proportion to the number of SNVs
    fig.write_image(
        str(output_dir / f"{patient_name}_DNA_heatmap.pdf"), 
        width = 500 + 10 * len(snv_df.index.tolist())
    )
    print(f'[INFO] Saved heatmap for {patient_name} to {output_dir / f"{patient_name}_DNA_heatmap.pdf"}.')
    fig.write_image(
        str(output_dir / f"{patient_name}_DNA_heatmap.png"), 
        width = 500 + 10 * len(snv_df.index.tolist()), scale=4,
    )
    print(f'[INFO] Saved heatmap for {patient_name} to {output_dir / f"{patient_name}_DNA_heatmap.png"}.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--patient_name', type=str, help='patient name')
    parser.add_argument('--patient_h5', type=str, help='input path to patient-wide h5')
    parser.add_argument('--snv_f', type=str, help='input path to SNV file')
    parser.add_argument('--clone_assignment', type=str, help='input path to clone assignment file')
    parser.add_argument('--output_dir', type=str, help='output directory')
    parser.add_argument('--post_condor', default=False, action='store_true', help='whether the clone assignment file is post-condor')
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    main(args)

