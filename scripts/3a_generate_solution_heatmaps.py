import os
import sys
import pandas as pd
import numpy as np
import math
import shutil
import argparse
import itertools

import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

def generate_condor_solution_heatmap(df_vaf, df_variant_readcounts, df_character_matrix, df_multistate, sorted_mutation_list, snv_annotations, germline_mutations, somatic_mutations):
    df_solution = df_variant_readcounts.copy()
    clustering = df_character_matrix['cluster_id']
    cell_order = list(clustering.sort_values().index)
    
    df_solution = df_solution.loc[cell_order]
    clustering = clustering[cell_order]

    for mutation in df_vaf.columns[:-1]:
        for cluster_id in df_multistate.index:
            df_solution.loc[clustering == cluster_id, mutation] = df_multistate.loc[cluster_id, mutation]
    
    #myColors = ('#636363',  '#334AFF', '#C1C1C1', '#FF3333')
    myColors = ('#FFFDAF', '#334AFF', '#C1C1C1', '#FF3333')
    cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))
    print(clustering)
    print('check')
    row_col_list = list(clustering.apply(lambda x: sns.color_palette()[int(x)]))
    
    df_solution = df_solution.drop(columns=['Unnamed: 0'])
    df_solution = df_solution[sorted_mutation_list]
    
    def mut_replace(x):
        x = x.replace(":", "_").replace("/" , "_").split('_')
        x[2], x[3] = x[3], x[2]
        return "_".join(x)
    
    df_solution.rename(columns=snv_annotations, inplace=True)
    _fig_h = 10
    _fig_w = len(sorted_mutation_list) * 0.5 + 2
    gs = sns.clustermap(df_solution, row_cluster=False, col_cluster=False, xticklabels=1,
                       cmap=cmap, figsize=(_fig_w,_fig_h), vmin=0, vmax=3, row_colors=row_col_list)

    gs.ax_heatmap.set_ylabel('cells', fontsize=20)
    gs.ax_heatmap.set_xlabel('mutations', fontsize=20)
    gs.ax_heatmap.axes.set_yticklabels([])
    gs.ax_heatmap.axes.set_xticklabels(gs.ax_heatmap.axes.get_xticklabels(), size = 10)

    # snv_annotations = {mut_replace(k): v for k, v in snv_annotations.items()}
    
    for tick_label in gs.ax_heatmap.axes.get_xticklabels():
        if tick_label.get_text() in germline_mutations:
            tick_label.set_color('green')
        elif tick_label.get_text() in somatic_mutations:
            tick_label.set_color('red')
        else: 
            tick_label.set_color('black')

    gs.ax_col_dendrogram.set_visible(False)
    gs.ax_row_dendrogram.set_visible(False)
    gs.ax_heatmap.tick_params(right=False)
    
    return gs
       
def generate_vaf_heatmap(df_vaf, df_character_matrix, df_variant_readcounts, df_total_readcounts, sorted_mutation_list, snv_annotations, germline_mutations, somatic_mutations): 

    clustering = df_character_matrix['cluster_id']
    cell_order = list(clustering.sort_values().index)

    df_vaf = df_vaf.loc[cell_order]
    df_variant_readcounts = df_variant_readcounts.loc[cell_order]
    df_total_readcounts = df_total_readcounts.loc[cell_order]
    clustering = clustering[cell_order]

    row_col_list = list(clustering.apply(lambda x: sns.color_palette()[int(x)]))
    
    df_vaf = df_vaf.drop(columns=['Unnamed: 0', 'cluster_id'])
    df_vaf = df_vaf[sorted_mutation_list]
    
    df_vaf = df_vaf.fillna(-1)
    
    #cmap = sns.color_palette("coolwarm", as_cmap=True)
    colors = ['#C1C1C1', '#334AFF', '#FF3333'] # first color is black, last is red
    cmap = LinearSegmentedColormap.from_list("Custom", colors, N=100)

    cmap.set_under('yellow')
    
    df_vaf.rename(columns=snv_annotations, inplace=True)
    
    print(df_vaf.columns)
    _fig_h = 10
    _fig_w = len(sorted_mutation_list) * 0.5 + 2
    gs = sns.clustermap(df_vaf, row_cluster=False, col_cluster=False, vmin=0, vmax=1, xticklabels=1,
                               cmap=cmap, figsize=(_fig_w,_fig_h), row_colors=row_col_list)
    
    gs.ax_heatmap.set_ylabel('cells', fontsize=20)
    gs.ax_heatmap.set_xlabel('mutations', fontsize=20)
    gs.ax_heatmap.axes.set_yticklabels([])
    gs.ax_heatmap.axes.set_xticklabels(gs.ax_heatmap.axes.get_xticklabels(), size = 10)

    def mut_replace(x):
        x = x.replace(":", "_").replace("/" , "_").split('_')
        x[2], x[3] = x[3], x[2]
        return "_".join(x)
    # snv_annotations = {mut_replace(k): v for k, v in snv_annotations.items()}
    print("==============")
    print(snv_annotations)
    print("==============")
    
    for tick_label in gs.ax_heatmap.axes.get_xticklabels():
        print(tick_label.get_text())
        if tick_label.get_text() in germline_mutations:
            tick_label.set_color('green')
        elif tick_label.get_text() in somatic_mutations:
            tick_label.set_color('red')
        else:
            tick_label.set_color('black')

    gs.ax_col_dendrogram.set_visible(False)
    gs.ax_row_dendrogram.set_visible(False)
    gs.ax_heatmap.tick_params(right=False)

    return gs

def main(args):
    df_character_matrix = pd.read_csv(args.c)
    df_multistate = pd.read_csv(args.s)

    df_vaf = pd.read_csv(args.v)
    df_variant_readcounts = pd.read_csv(args.a)
    df_total_readcounts = pd.read_csv(args.t)
    df_amplicon = pd.read_csv(args.i)
    if not "chr" not in df_amplicon.columns:
        if not "chrom" in df_amplicon.columns:
            raise IndexError("amplicon metadata file must contain the column named 'chrom' or 'chr'")
        else:
            # rename
            df_amplicon = df_amplicon.rename(columns={"chrom": "chr"})

    selected_mutation_list = [col for col in df_vaf.columns if col.startswith('chr')]
    
    mut_data = []
    for mutation in selected_mutation_list:
        mut_chr = mutation.split(':')[0].lstrip('chr')
        print(mutation)
        mut_pos = int(mutation.split(':')[1]) # <<<<<<<<<<
        mut_gene = df_amplicon[((df_amplicon['chrom'] == mut_chr) &
                                (df_amplicon['min_pos'] <= mut_pos) &
                                (df_amplicon['max_pos'] >= mut_pos-1))]['gene'].values[0]
        mut_data.append([mutation, mut_chr, mut_pos, mut_gene])
    
    df_mut_meta_selected = pd.DataFrame(mut_data, columns = ['mutation', 'chr', 'pos', 'gene'])
    df_mut_meta_selected['chr'] = df_mut_meta_selected['chr'].apply(lambda x: int(x) if x != 'X' else 23)
    df_mut_meta_selected = df_mut_meta_selected.sort_values(['chr', 'pos'])
    sorted_mutation_list = list(df_mut_meta_selected['mutation'])
    

    # @HZ use annotated SNV file instead of CRAVAT
    snvs = pd.read_csv(args.snvs, sep="\t", comment="#")
    snvs.replace(np.nan, "", inplace=True)
    snvs.set_index("condensed_format", inplace=True)
    snv_annotations = snvs["HGVSp"].to_dict()

    germline_mutations_fp = args.g
    germline_mutations = []
    with open(germline_mutations_fp, "r") as file:
        for line in file:
            # Remove leading/trailing whitespace and newline characters
            mut = line.strip()
            germline_mutations.append(mut)
    germline_mutations = [snv_annotations[mut] for mut in germline_mutations]

    somatic_mutations_fp = args.m
    somatic_mutations = []
    with open(somatic_mutations_fp, "r") as file:
        for line in file:
            # Remove leading/trailing whitespace and newline characters
            mut = line.strip()
            somatic_mutations.append(mut)
    somatic_mutations = [snv_annotations[mut] for mut in somatic_mutations]

    gs_sol = generate_condor_solution_heatmap(df_vaf, df_variant_readcounts, df_character_matrix, df_multistate, sorted_mutation_list, snv_annotations, germline_mutations, somatic_mutations)
    gs_vaf = generate_vaf_heatmap(df_vaf, df_character_matrix, df_variant_readcounts, df_total_readcounts, sorted_mutation_list, snv_annotations, germline_mutations, somatic_mutations)

    gs_sol.savefig(args.o) 
    gs_vaf.savefig(args.p) 

    return 0

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=str, help='dataset name')
    parser.add_argument('-c', type=str, help='input path to character matrix')
    parser.add_argument('-v', type=str, help='input path to VAF matrix')
    parser.add_argument('-s', type=str, help='input path to CONDOR solution matrix')
    parser.add_argument('-g', type=str, help='input path to CONDOR input germline mutations list')
    parser.add_argument('-m', type=str, help='input path for somatic mutation list')
    parser.add_argument('-a', type=str, help='input path to alternate readcount matrix')
    parser.add_argument('-t', type=str, help='input path to total readcount matrix')
    parser.add_argument('-i', type=str, help='input path to amplicon panel')
    # parser.add_argument('-l', type=str, help='input path to CRAVAT files directory')
    parser.add_argument('-o', type=str, help='output path to CONDOR solution heatmap png')
    parser.add_argument('-p', type=str, help='output path to VAF heatmap png')
    parser.add_argument(
        "-snvs", type=str, default=None, help="input path to manually annotated SNVs"
    )
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    main(args)
