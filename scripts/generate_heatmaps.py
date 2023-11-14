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

def generate_condor_solution_heatmap(df_vaf, df_variant_readcounts, df_character_matrix, df_multistate, sorted_mutation_list, cravat_df, germline_mutations, somatic_mutations):
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


    mapping_dict = list(map(mut_replace, cravat_df['var_ann'].index))

    mapping_dict = dict(zip(cravat_df['var_ann'].index, mapping_dict))

    # Rename the index using the new_index
    cravat_df.rename(index=mapping_dict, inplace=True)
    
    mapping_dict = {idx:cravat_df['var_ann'].loc[idx].values[0] for idx in cravat_df['var_ann'].index}
    
    df_solution.rename(columns=mapping_dict, inplace=True)

    gs = sns.clustermap(df_solution, row_cluster=False, col_cluster=False, xticklabels=1,
                       cmap=cmap, figsize=(20,10), vmin=0, vmax=3, row_colors=row_col_list)

    gs.ax_heatmap.set_ylabel('cells', fontsize=20)
    gs.ax_heatmap.set_xlabel('mutations', fontsize=20)
    gs.ax_heatmap.axes.set_yticklabels([])
    gs.ax_heatmap.axes.set_xticklabels(gs.ax_heatmap.axes.get_xticklabels(), size = 10)
    

    inv_mapping_dict = {v: k for k, v in mapping_dict.items()}
    
    for tick_label in gs.ax_heatmap.axes.get_xticklabels():
        if inv_mapping_dict[tick_label.get_text()] in germline_mutations:
            tick_text = tick_label.get_text()
            tick_label.set_color('green')
        elif inv_mapping_dict[tick_label.get_text()] in somatic_mutations:
            tick_text = tick_label.get_text()
            tick_label.set_color('red')
        else: 
            tick_text = tick_label.get_text()
            tick_label.set_color('black')

    gs.ax_col_dendrogram.set_visible(False)
    gs.ax_row_dendrogram.set_visible(False)
    gs.ax_heatmap.tick_params(right=False)
    
    return gs
    

    
def generate_vaf_heatmap(df_vaf, df_character_matrix, df_variant_readcounts, df_total_readcounts, sorted_mutation_list, cravat_df, germline_mutations, somatic_mutations): 

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
    
    def mut_replace(x):
        x = x.replace(":", "_").replace("/" , "_").split('_')
        #x[2], x[3] = x[3], x[2]
        return "_".join(x)


    mapping_dict = list(map(mut_replace, cravat_df['var_ann'].index))

    mapping_dict = dict(zip(cravat_df['var_ann'].index, mapping_dict))

    # Rename the index using the new_index
    cravat_df.rename(index=mapping_dict, inplace=True)
    
    mapping_dict = {idx:cravat_df['var_ann'].loc[idx].values[0] for idx in cravat_df['var_ann'].index}
    
    df_vaf.rename(columns=mapping_dict, inplace=True)
    
    print(df_vaf.columns)
    gs = sns.clustermap(df_vaf, row_cluster=False, col_cluster=False, vmin=0, vmax=1, xticklabels=1,
                               cmap=cmap, figsize=(20,10), row_colors=row_col_list)
    
    gs.ax_heatmap.set_ylabel('cells', fontsize=20)
    gs.ax_heatmap.set_xlabel('mutations', fontsize=20)
    gs.ax_heatmap.axes.set_yticklabels([])
    gs.ax_heatmap.axes.set_xticklabels(gs.ax_heatmap.axes.get_xticklabels(), size = 10)

    inv_mapping_dict = {v: k for k, v in mapping_dict.items()}
    
    for tick_label in gs.ax_heatmap.axes.get_xticklabels():
        if inv_mapping_dict[tick_label.get_text()] in germline_mutations:
            tick_text = tick_label.get_text()
            tick_label.set_color('green')
        elif inv_mapping_dict[tick_label.get_text()] in somatic_mutations:
            tick_label.set_color('red')
        else:
            tick_text = tick_label.get_text()
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
        if "chrom" in df_amplicon.columns:
            raise IndexError("amplicon metadata file must contain the column named 'chrom' or 'chr'")
        else:
            # rename
            df_amplicon = df_amplicon.rename(columns={"chrom": "chr"})
    local_directory = args.l
    dataset = args.d


    selected_mutation_list = [col for col in df_vaf.columns if col.startswith('chr')]
    
    mut_data = []
    for mutation in selected_mutation_list:
        mut_chr = mutation.split('_')[0].lstrip('chr')
        mut_pos = int(mutation.split('_')[1])
        mut_gene = df_amplicon[((df_amplicon['chrom'] == mut_chr) &
                                (df_amplicon['min_pos'] <= mut_pos) &
                                (df_amplicon['max_pos'] >= mut_pos-1))]['gene'].values[0]
        mut_data.append([mutation, mut_chr, mut_pos, mut_gene])
    
    df_mut_meta_selected = pd.DataFrame(mut_data, columns = ['mutation', 'chr', 'pos', 'gene'])
    df_mut_meta_selected['chr'] = df_mut_meta_selected['chr'].apply(lambda x: int(x) if x != 'X' else 23)
    df_mut_meta_selected = df_mut_meta_selected.sort_values(['chr', 'pos'])
    sorted_mutation_list = list(df_mut_meta_selected['mutation'])
    

    cravat_df = None
    for file in os.listdir(local_directory + dataset + '/'):
        if file.startswith(dataset):
        # CRAVAT File
            if file.endswith('.txt'):
                cravat_f = local_directory + dataset + '/' + file
                cravat_df = pd.read_csv(cravat_f, sep='\t', index_col=0, header=[0,1])
    
    germline_mutations_fp = args.g
    germline_mutations = []
    
    with open(germline_mutations_fp, "r") as file:
        for line in file:
            # Remove leading/trailing whitespace and newline characters
            mut = line.strip()
            germline_mutations.append(mut)
    
    somatic_mutations_fp = args.m
    somatic_mutations = []
    
    with open(somatic_mutations_fp, "r") as file:
        for line in file:
            # Remove leading/trailing whitespace and newline characters
            mut = line.strip()
            somatic_mutations.append(mut)

    gs_sol = generate_condor_solution_heatmap(df_vaf, df_variant_readcounts, df_character_matrix, df_multistate, sorted_mutation_list, cravat_df, germline_mutations, somatic_mutations)
    gs_vaf = generate_vaf_heatmap(df_vaf, df_character_matrix, df_variant_readcounts, df_total_readcounts, sorted_mutation_list, cravat_df, germline_mutations, somatic_mutations)

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
    parser.add_argument('-l', type=str, help='input path to CRAVAT files directory')
    parser.add_argument('-o', type=str, help='output path to CONDOR solution heatmap png')
    parser.add_argument('-p', type=str, help='output path to VAF heatmap png')
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    main(args)
