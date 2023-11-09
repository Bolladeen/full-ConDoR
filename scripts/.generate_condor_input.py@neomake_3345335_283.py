import os
import sys
import pandas as pd
import numpy as np
import math
import shutil
import argparse
import itertools

import h5py
import mosaic.io as mio
import plotly.express as px
from tea.cravat import get_technical_artifact_mask
from tea.format import isNaN
from tea.plots import plot_snv_clone


def get_filtered_mutations(cravat_df, sample_obj, sample_name):
    num_cells = sample_obj.dna.shape[0]

    # all params are set to default
    technical_mask = get_technical_artifact_mask(
        cravat_df,
        num_cells = num_cells,
        bq_prev_threshold = None, #0.005,
        normals_pon_occurence = 4, #4
        rescue_1000genome_af = 0.01,
        filter_broad_wes_pon = False
        )

    #technical_mask = technical_mask & (cravat_df[('Tapestri_result', 'sc_mut_prev')] >= 0.01 * num_cells)
    technical_mask =  (cravat_df[('bulk_comparison', 'bulk-matched_bulk_normal-AF')] > 0) | (cravat_df[('bulk_comparison', 'bulk-matched_bulk_cohort-AF')] > 0)

    bulk_germline_snv_set = set(cravat_df[(cravat_df[('bulk_comparison', 'bulk-matched_bulk_normal-AF')] > 0) & technical_mask].index)

    bulk_somatic_snv_set = set(cravat_df[(cravat_df[('bulk_comparison', 'bulk-matched_bulk_cohort-AF')] > 0) & technical_mask].index)
    voi = cravat_df[technical_mask].index.tolist()

    ann = cravat_df.loc[voi, :].index.map(
    lambda x:
        cravat_df.loc[x, ('Variant Annotation', 'Gene')] + ' ' + cravat_df.loc[x, ('Variant Annotation', 'Protein Change')] if not isNaN(cravat_df.loc[x, ('Variant Annotation','Protein Change')])
        else cravat_df.loc[x, ('Variant Annotation','Gene')] + ' ' + cravat_df.loc[x, ('Variant Annotation','Sequence Ontology')]
    )

    ann_map = dict(zip(voi, ann))

    total_mutations = []
    germline_mutations = []
    for var_i in ann_map:
        total_mutations.append(var_i)
        # germline
        if var_i in bulk_germline_snv_set:
            germline_mutations.append(var_i)
        else:
            pass

    return total_mutations, germline_mutations

def generate_condor_input(sample_name, cn_assignment_df, args, bin_thres=0.5):
    merged_cn_assignment_df = cn_assignment_df.copy()
    merged_cn_assignment_df = merged_cn_assignment_df.rename(columns={'Unnamed: 0':'sample'})
    merged_cn_assignment_df['cell_barcode'] =  merged_cn_assignment_df['sample'].astype(str) + ':' + merged_cn_assignment_df['cell_barcode'].astype(str)
    merged_cn_assignment_df.drop(['sample'], axis=1, inplace=True)

    df_total_samples = []
    df_alt_samples = []
    common_mutations = []
    germline_mutations_list = []
    somatic_mutations = []
    cravat_f = args.l + sample_name + '/' + sample_name + '_CRAVAT_output_cleaned.txt'
    cravat_df = pd.read_csv(cravat_f, sep='\t', index_col=0, header=[0,1])
    for file in os.listdir(args.l + sample_name):
        if file.startswith(sample_name):
            if file.endswith('.h5'):
                hdf5_f = args.l + sample_name + '/' + file
                sample = file.split('.')[0]
                sample_obj = mio.load(hdf5_f)
                
                df_alt_snv = sample_obj.dna.get_attribute('alt_read_count', constraint='row')
                df_total_snv = sample_obj.dna.get_attribute('DP', constraint='row')
                print(df_total_snv.shape)
                cmuts, gmuts = get_filtered_mutations(cravat_df, sample_obj, sample_name)
                smuts = set(cmuts) - set(gmuts)
                common_mutations.append(set(cmuts))
                germline_mutations_list.append(set(gmuts))
                somatic_mutations.append(smuts)

                df_alt_snv.reset_index(inplace=True)
                df_alt_snv = df_alt_snv.rename(columns = {'index':'cell_barcode'})

                df_alt_samples.append(df_alt_snv)
                    
                df_total_snv.reset_index(inplace=True)
                df_total_snv = df_total_snv.rename(columns = {'index':'cell_barcode'})
                
                df_total_samples.append(df_total_snv)
   
    germline_mutations_list = set.intersection(*germline_mutations_list)
    somatic_mutations = set.union(*somatic_mutations)
    # RA17_13
    #somatic_mutations = set(['chr9:21974695:G/GT', 'chr12:25398285:C/A', 'chr18:48581173:G/T', 'chr9:101891277:C/T', 'chr9:101900209:C/T', 'chr17:7579717:GA/G'])

    common_mutations = list(set(list(somatic_mutations) + list(germline_mutations_list)))

    def mut_replace(x):
        x = x.replace(":", "_").replace("/" , "_").split('_')
        x[2], x[3] = x[3], x[2]
        return "_".join(x)

    common_mutations = list(map(mut_replace, common_mutations))
    germline_mutations_list = list(map(mut_replace, germline_mutations_list))
    print(len(common_mutations))
    
    df_total = pd.concat(df_total_samples, ignore_index=True)
    df_total = df_total.rename(columns = {'DP':'cell_barcode'})
    df_total = df_total.rename(columns={c: mut_replace(c) for c in df_total.columns if c not in ['cell_barcode']})

    df_total = df_total[['cell_barcode'] + common_mutations]
    df_total = df_total.set_index('cell_barcode')
    df_total = df_total.fillna(0)
    
    df_alt = pd.concat(df_alt_samples)
    df_alt = df_alt.rename(columns = {'alt_read_count':'cell_barcode'})
    df_alt = df_alt.rename(columns={c: mut_replace(c) for c in df_alt.columns if c not in ['cell_barcode']})

    df_alt = df_alt[['cell_barcode'] + common_mutations]
    df_alt = df_alt.set_index('cell_barcode')
    df_alt = df_alt.fillna(0)
    
    print(df_total)
    print(df_alt)

    df_character_mat = df_total.copy()
    df_character_mat =  df_alt.divide(df_total)


    print(df_character_mat.index)
    print(merged_cn_assignment_df['cell_barcode'])
    
    def rename_barcode(s):
        return s.split(':')[1] + '-' + s.split(':')[0]
    
    merged_cn_assignment_df['cell_barcode'] = merged_cn_assignment_df['cell_barcode'].apply(rename_barcode)

    df_character_mat = pd.merge(df_character_mat, merged_cn_assignment_df,left_on=df_character_mat.index, right_on='cell_barcode', how='left')
    df_character_mat = df_character_mat.set_index('cell_barcode')
    
    df_character_mat.rename(columns={'clone_id': 'cluster_id'}, inplace=True)
    
    l_ids = list(df_character_mat['cluster_id'].unique())
    l_ids.sort()
    l_ids_dict = {}
    index = 0
    for i in l_ids:
        l_ids_dict[i] = index
        index += 1

    df_character_mat['cluster_id'] =  df_character_mat['cluster_id'].replace(l_ids_dict)

    df_character_vaf = df_character_mat.copy()
    df_character_mat[df_character_mat.columns[:-1]] = df_character_mat[df_character_mat.columns[:-1]].applymap(lambda x: -1 if pd.isna(x) else 0 if x < bin_thres else 1)

    return df_total, df_alt, df_character_mat, df_character_vaf, germline_mutations_list

def main(args):
    dataset = args.d
    cn_assignment_df = pd.read_csv(args.i)
    df_total, df_alt, df_character_mat, df_character_vaf, germline_mutations = generate_condor_input(dataset, cn_assignment_df, args)
    with open(args.g, 'w') as fp:
        for item in germline_mutations:
            fp.write("%s\n" % item)

    df_character_vaf.to_csv(args.v, index_label='')
    df_character_mat.to_csv(args.m, index_label='')
    df_alt.to_csv(args.a,index_label='')
    df_total.to_csv(args.t,index_label='')

    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=str, help='dataset name')
    parser.add_argument('-l', type=str, help='input path for hdf5 and CRAVAT files')
    parser.add_argument('-i', type=str, help='input path to refined cell assignments csv file')
    parser.add_argument('-v', type=str, help='output path for ConDoR VAF matrix')
    parser.add_argument('-m', type=str, help='output path for ConDoR binary matrix')
    parser.add_argument('-a', type=str, help='output path for ConDoR alternate readcount matrix')
    parser.add_argument('-t', type=str, help='output path for ConDoR total readcount matrix ')
    parser.add_argument('-g', type=str, help='output path for ConDoR germline mutation list')
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    main(args)
