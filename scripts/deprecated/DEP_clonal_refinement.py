import os
import sys
import pandas as pd
import numpy as np
import math
import shutil
import argparse
import itertools
import yaml

import h5py
import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter
from collections import defaultdict
import mosaic.io as mio
import plotly.express as px
from tea.cravat import get_technical_artifact_mask
from tea.format import isNaN
from tea.plots import plot_snv_clone
from pathlib import Path

def get_figure_and_germline_mutations(cravat_df, sample_obj, analysis_config, manually_annotated_snvs):
    num_cells = sample_obj.dna.shape[0]
    
    snv_selection_params = analysis_config['snv_selection_params']
        
    # 1. mutational prevalence
    mut_prev_threshold = snv_selection_params['mut_prev_threshold']
    if not type(mut_prev_threshold) is list:
        mut_prev_threshold = [mut_prev_threshold]

    # 2. technical artifact filters
    bq_prev_threshold = snv_selection_params['bq_prev_threshold']
    if bq_prev_threshold is not None and type(bq_prev_threshold) is not float: # single-value
        raise ValueError(f"bq_prev_threshold must be float, not {type(bq_prev_threshold)}")

    normals_occurences = snv_selection_params['normals_occurences'] if 'normals_occurences' in snv_selection_params else 3 # <--- default to 3

    if 'ado_threshold' in snv_selection_params and snv_selection_params['ado_threshold'] is not None:
        print(f"[INFO] filtering out SNVs with ADO > {snv_selection_params['ado_threshold']} in ANY sample...")
        ado_threshold = snv_selection_params['ado_threshold']
    else:
        ado_threshold = None

    # 3. functional SNVs
    topic = snv_selection_params['topic']
    if not type(topic) is str: # single-value
        raise ValueError(f"topic must be str, not {type(topic)}")
    try: 
        func_only = snv_selection_params['func_only']
        func_only = bool(func_only)
    except KeyError:
        func_only = False
    except TypeError:
        func_only = False

    # 4. germline snvs
    germline_attrs = {}
    for germline_attr_i in ['remove_hom_germline_snps_af', 'rescue_1000genome_af']:
        if not germline_attr_i in snv_selection_params:
            germline_attrs[germline_attr_i] = False
        else:
            germline_attrs[germline_attr_i] = snv_selection_params[germline_attr_i]

    
    if 'filter_TtoC_artifact' in snv_selection_params:
        try:
            filter_TtoC_artifact = snv_selection_params['filter_TtoC_artifact']['filter'] 
            filter_TtoC_artifact_lower_thres = snv_selection_params['filter_TtoC_artifact']['lower_thres']
            filter_TtoC_artifact_upper_thres = snv_selection_params['filter_TtoC_artifact']['upper_thres']
        except KeyError:
            raise ValueError(f"[ERROR] filter_TtoC_artifact must have keys ['filter', 'lower_thres', 'upper_thres']")
    else:
        if filter_TtoC_artifact_lower_thres >= filter_TtoC_artifact_upper_thres:
            raise ValueError(f"[ERROR] filter_TtoC_artifact_lower_thres must be strictly smaller than filter_TtoC_artifact_upper_thres")
        else:
            if mut_prev_i < 0.01:
                print(f"[INFO] mut_prev_i is lower than default upper_thres (0.01) for T>C filter. The filter will be applied.")
                filter_TtoC_artifact = True
                filter_TtoC_artifact_lower_thres = mut_prev_threshold[0]
                filter_TtoC_artifact_upper_thres = 0.01 
            else:
                print(f"[WARNING] mut_prev_i is higher than default upper_thres (0.01) for T>C filter. The filter will not be applied.")
                filter_TtoC_artifact = False

    voi_union = set()
    voi_count_union = {}
    ann_map_union = {}
    bulk_germline_vars = set()
    bulk_somatic_vars = set()
    TtoC_artifact_blacklist = set()
    ado_blacklist = set()
    
    mask = get_technical_artifact_mask(
            cravat_df,
            num_cells = num_cells, 
            bq_prev_threshold = bq_prev_threshold,
            normals_pon_occurence=normals_occurences,
            rescue_1000genome_af = germline_attrs['rescue_1000genome_af'],
            filter_broad_wes_pon = False)
    
    if filter_TtoC_artifact:
        tc_mask = (cravat_df.index.str.endswith('T/C')) & (cravat_df[('Tapestri_result', 'sc_mut_prev')] >= filter_TtoC_artifact_lower_thres * num_cells) & (cravat_df[('Tapestri_result', 'sc_mut_prev')] <= filter_TtoC_artifact_upper_thres * num_cells)
        mask = mask & ~tc_mask
        TtoC_artifact_blacklist = TtoC_artifact_blacklist.union(
        set(cravat_df.index[tc_mask].tolist())
        )

    # filters on functional SNVs
    if func_only:
        mask = mask & ~cravat_df[('Variant Annotation', 'Sequence Ontology')].isin(NONFUNC_SO)

    voi = cravat_df.index[mask].tolist()

    # filters on mut_prev_threshold
    prev_filtered_vars = sample_obj.dna.ids()[
        sample_obj.dna.get_attribute("mut_filtered", constraint="row").sum(axis=0) >= (mut_prev_threshold[0] * num_cells)
    ]
    # take intersection
    voi = [ v for v in voi if v in prev_filtered_vars ]
    

    if ado_threshold is not None:
        ado_high_vars = sample_obj.dna.ids()[
            (sample_obj.dna.get_attribute('NGT',constraint='row') == 3).sum(axis=0) > (ado_threshold*num_cells)
        ]
        voi = [ v for v in voi if v not in ado_high_vars ]
        ado_blacklist = ado_blacklist.union(set(ado_high_vars))

    voi_mut_prev = Counter(cravat_df.loc[voi, ('Tapestri_result', 'sc_mut_prev')].to_dict())

    ann = cravat_df.loc[voi, :].index.map(
        lambda x: 
        cravat_df.loc[x, ('Variant Annotation', 'Gene')] + ' ' + cravat_df.loc[x, ('Variant Annotation', 'Protein Change')] if not isNaN(cravat_df.loc[x, ('Variant Annotation','Protein Change')])
        else cravat_df.loc[x, ('Variant Annotation','Gene')] + ' ' + cravat_df.loc[x, ('Variant Annotation','Sequence Ontology')]
    )
    ann_map = dict(zip(voi, ann))


    print(len(voi))

    try:
        bulk_normal_vars = cravat_df.index[(cravat_df[('bulk_comparison', 'bulk-matched_bulk_normal-AF')] > 0)]
    except KeyError:
        print(f'bulk normal annotation not found in CRAVAT DF')
    else:
        if len(bulk_normal_vars) == 0:
            print(f'[WARNING] No bulk normal SNVs detected')
        bulk_germline_vars.update(bulk_normal_vars)

    # select SNVs detected in matched bulk cohort
    try:
        bulk_cohort_vars = cravat_df.index[(cravat_df[('bulk_comparison', 'bulk-matched_bulk_cohort-AF')] > 0)]
    except KeyError:
        print(f'bulk tumor annotation not found in CRAVAT DF')
    else:
        if len(bulk_cohort_vars) == 0:
            print(f'[WARNING] No bulk cohort SNVs detected')
        bulk_somatic_vars.update(bulk_cohort_vars)
    

    voi_union = voi_union.union(set(voi))
    voi_count_union.update(voi_mut_prev)
    ann_map_union.update(ann_map)

    # remove SNVs that are blacklisted
    print(f"[DEBUG] {len(voi_union)} SNVs before blacklist filtering")
    voi_union = voi_union.difference(TtoC_artifact_blacklist)
    print(f"[DEBUG] {len(voi_union)} SNVs after TtoC blacklist filtering")
    voi_union = voi_union.difference(ado_blacklist)
    print(f"[DEBUG] {len(voi_union)} SNVs after ADO blacklist filtering")
    
    total_mutations = []
    germline_mutations = []
    for var_i in ann_map:
        total_mutations.append(var_i)
        # germline
        overall_af = sample_obj.dna.get_attribute('alt_read_count', constraint='row').sum(axis=0) / sample_obj.dna.get_attribute('DP', constraint='row').sum(axis=0)
        #if overall_af[var_i] >= snv_selection_params['remove_hom_germline_snps_af']:
            #germline_mutations.append(var_i)
        if var_i in bulk_germline_vars:
            germline_mutations.append(var_i)
        else:
            pass
    
    fig = None
    
    if manually_annotated_snvs is not None and os.stat(manually_annotated_snvs).st_size != 0:
        snvs = pd.read_csv(manually_annotated_snvs, sep='\t', comment='#')
        snvs.replace(np.nan, '', inplace=True)
        snvs = snvs[snvs['annotation'].str.startswith('germline')]['condensed_format'].to_list()
        germline_mutations = snvs

    germline_mutations = list(set(germline_mutations))
    print(germline_mutations)
    return fig, germline_mutations


def binarize_NGT_matrix(NGT_df):
    NGT_df[NGT_df.columns[:-2]] = np.where(NGT_df[NGT_df.columns[:-2]] == 1, 0, NGT_df[NGT_df.columns[:-2]])
    NGT_df[NGT_df.columns[:-2]] = np.where(NGT_df[NGT_df.columns[:-2]] == 2, 1, NGT_df[NGT_df.columns[:-2]])

    return NGT_df

def profile_distance(medians, cell, cell_barcode, cell_profile):
    distances = {}
    for median_clone, median_profile in medians.items():
        mask = median_profile != 3
        mask2 = cell_profile != 3
        mask = np.logical_and(mask, mask2)
        series1 = median_profile[mask]
        series2 = cell_profile[mask]
        absolute_sum = np.sqrt(((series1 - series2)**2).sum())
        distances[median_clone] = absolute_sum/series2.shape[0]

    return distances


def median_distance(medians, clone1, clone2):
    median_profile1 = medians[clone1]
    median_profile2 = medians[clone2]
    mask = median_profile1 != 3
    mask2 = median_profile2 != 3
    mask = np.logical_and(mask, mask2)
    series1 = median_profile1[mask]
    series2 = median_profile2[mask]
    absolute_sum = np.sqrt(((series1 - series2)**2).sum())
    distance = absolute_sum/series2.shape[0]

    return distance


def generate_median_cluster_profiles(NGT_df):
    median_clusters = {}
    median_clusters_all = {}

    for clone in NGT_df['clone_id'].unique():
        NGT_curr_clone = NGT_df[NGT_df['clone_id'] == clone][NGT_df.columns[:-2]]
        if NGT_curr_clone.shape[0] > 0:
            
            if NGT_curr_clone.mode().squeeze(axis=0).shape[0] == 2 or NGT_curr_clone.mode().squeeze(axis=0).shape[0] == 3:
                median_clusters_all[clone] = NGT_curr_clone.mode().squeeze(axis=0).iloc[0]
            else:
                median_clusters_all[clone] = NGT_curr_clone.mode().squeeze(axis=0)

        if NGT_curr_clone.shape[0] > 0.05 * NGT_df.shape[0]:
            if NGT_curr_clone.mode().squeeze(axis=0).shape[0] == 2 or NGT_curr_clone.mode().squeeze(axis=0).shape[0] == 3:
                median_clusters[clone] = NGT_curr_clone.mode().squeeze(axis=0).iloc[0]
            else:
                median_clusters[clone] = NGT_curr_clone.mode().squeeze(axis=0)

    return median_clusters, median_clusters_all


def copy_number_distance(df_tapestri_clones, clone1, clone2):
    c1_profile = np.array(df_tapestri_clones.loc[clone1].tolist())
    c2_profile = np.array(df_tapestri_clones.loc[clone2].tolist())
    absolute_sum = np.abs(c1_profile - c2_profile).sum()
    absolute_sum = absolute_sum/len(c2_profile)
    return absolute_sum

def combine_similar_final_clusters(NGT_df, cn_assignment_df, median_clusters, df_tapestri_clones):
    normal_idx = -1
    if len(df_tapestri_clones.index[df_tapestri_clones.eq(2.0).all(axis=1)].to_list()) > 0:
        normal_idx = df_tapestri_clones.index[df_tapestri_clones.eq(2.0).all(axis=1)].to_list()[0]
    for c1 in NGT_df['clone_id'].unique():
        for c2 in NGT_df['clone_id'].unique():
            if c1 != c2:
                if median_distance(median_clusters, c1, c2) <= 0.10:
                    if copy_number_distance(df_tapestri_clones, c1, c2) <= 0.40:
                        c1_size = NGT_curr_clone = NGT_df[NGT_df['clone_id'] == c1][NGT_df.columns[:-1]].shape[0]
                        c2_size = NGT_curr_clone = NGT_df[NGT_df['clone_id'] == c2][NGT_df.columns[:-1]].shape[0]
                            
                        if c1 != normal_idx and c2 != normal_idx:
                            if c1_size >= c2_size:
                                cn_assignment_df['clone_id'] = cn_assignment_df['clone_id'].replace(c2, c1)
                                NGT_df['clone_id'] = NGT_df['clone_id'].replace(c2, c1)
                            else:
                                cn_assignment_df['clone_id'] = cn_assignment_df['clone_id'].replace(c1, c2)
                                NGT_df['clone_id'] = NGT_df['clone_id'].replace(c1, c2)

                        else:
                            if c1 == normal_idx:
                                cn_assignment_df['clone_id'] = cn_assignment_df['clone_id'].replace(c2, c1)
                                NGT_df['clone_id'] = NGT_df['clone_id'].replace(c2, c1)
                            else:
                                cn_assignment_df['clone_id'] = cn_assignment_df['clone_id'].replace(c1, c2)
                                NGT_df['clone_id'] = NGT_df['clone_id'].replace(c1, c2)


    return NGT_df, cn_assignment_df, df_tapestri_clones

def merge_large_clusters(NGT_df, NGT_df_homvaf, cn_assignment_df, median_clusters,  count, sample_name):
    for clone in NGT_df['clone_id'].unique():
        NGT_curr_clone = NGT_df[NGT_df['clone_id'] == clone][NGT_df.columns[:-1]]
        if NGT_curr_clone.shape[0] > 0.05 * NGT_df.shape[0]:
            for cell, cell_profile in NGT_curr_clone.iterrows():
                cell_barcode = cell_profile['cell_barcode']
                cell_profile = cell_profile.drop('cell_barcode')
                distances = profile_distance(median_clusters, cell, cell_barcode, cell_profile)
                distances = {k:v/distances[clone] for k,v in distances.items()}
                k,v = min(distances.items(), key=lambda x: x[1])
                if v < 0.85:
                    count += 1
                    cn_assignment_df.loc[(cn_assignment_df.index == cell) & (cn_assignment_df['cell_barcode'] == cell_barcode), 'clone_id'] = k
                    NGT_df_homvaf.loc[(NGT_df_homvaf.index == cell) & (NGT_df_homvaf['cell_barcode'] == cell_barcode), 'clone_id'] = k
                    NGT_df.loc[(NGT_df.index == cell) & (NGT_df['cell_barcode'] == cell_barcode), 'clone_id'] = k

    return NGT_df, cn_assignment_df, count

def dissolve_small_clusters(NGT_df, NGT_df_homvaf, cn_assignment_df, median_clusters,  count, sample_name):
    for clone in NGT_df['clone_id'].unique():
        NGT_curr_clone = NGT_df[NGT_df['clone_id'] == clone][NGT_df.columns[:-1]]
        if NGT_curr_clone.shape[0] > 0 and NGT_curr_clone.shape[0] <= 0.05 * NGT_df.shape[0]:
            for cell, cell_profile in NGT_curr_clone.iterrows():
                cell_barcode = cell_profile['cell_barcode']
                cell_profile = cell_profile.drop('cell_barcode')
                distances = profile_distance(median_clusters, cell, cell_barcode, cell_profile)
                distances = {k:v for k,v in distances.items()}
                distances.pop(clone, None)
                k,v = min(distances.items(), key=lambda x: x[1])
                count += 1
                cn_assignment_df.loc[(cn_assignment_df.index == cell) & (cn_assignment_df['cell_barcode'] == cell_barcode), 'clone_id'] = k
                NGT_df_homvaf.loc[(NGT_df_homvaf.index == cell) & (NGT_df_homvaf['cell_barcode'] == cell_barcode), 'clone_id'] = k
                NGT_df.loc[(NGT_df.index == cell) & (NGT_df['cell_barcode'] == cell_barcode), 'clone_id'] = k

    return NGT_df, cn_assignment_df, count




def split_sample_clones(sample_name, NGT_df, NGT_df_homvaf, df_tapestri_clones, cn_assignment_df, added_clone_mut_dict):
    count = 0
    added_clone = False
    
    median_clusters, median_clusters_all = generate_median_cluster_profiles(NGT_df)
    NGT_df, cn_assignment_df, count = merge_large_clusters(NGT_df, NGT_df_homvaf, cn_assignment_df, median_clusters,  count, sample_name)
    NGT_df, cn_assignment_df, count = dissolve_small_clusters(NGT_df, NGT_df_homvaf,  cn_assignment_df, median_clusters,  count, sample_name)
    #split

    max_cn = cn_assignment_df['clone_id'].max()
    for clone in NGT_df['clone_id'].unique():
        candidates_not_found = False
        while(candidates_not_found == False):
            candidate_muts = {}
            NGT_curr_clone = NGT_df_homvaf[NGT_df_homvaf['clone_id'] == clone][NGT_df_homvaf.columns[:-1]]
            for mut in NGT_curr_clone.columns:
                if mut != 'cell_barcode':
                    non_missing_NGT_curr_clone = NGT_curr_clone[mut].replace(3, np.NaN)
                    mean = non_missing_NGT_curr_clone.mean()
                    if mean < 0.75:
                        if 2 in NGT_curr_clone[mut].unique():
                            if NGT_curr_clone[mut].value_counts()[2.0] > 0.10 * NGT_df_homvaf.shape[0] and NGT_curr_clone[mut].value_counts()[3.0] < 0.30 * NGT_curr_clone.shape[0]:
                                candidate_muts[mut] = NGT_curr_clone[mut].value_counts()[2.0]/NGT_curr_clone.shape[0]

            candidate_muts = {k: v for k, v in sorted(candidate_muts.items(), key=lambda item: item[1], reverse=True)}
            if len(candidate_muts) < 2:
                candidates_not_found = 1
            else:
                added_clone = True
                max_cn = max_cn + 1
                candidate = list(candidate_muts.keys())[0]
                print(candidate, clone, sample_name)
                added_clone_mut_dict[max_cn] = candidate
                new_clone_profile = df_tapestri_clones.loc[clone].copy()
                df_tapestri_clones = pd.concat([df_tapestri_clones, pd.DataFrame([new_clone_profile])])
                #df_tapestri_clones = df_tapestri_clones.append(new_clone_profile)
                df_tapestri_clones.index.values[-1] = max_cn
                for cell, cell_profile in NGT_curr_clone.iterrows():
                    cell_barcode = cell_profile['cell_barcode']
                    cell_profile = cell_profile.drop('cell_barcode')
                    if cell_profile.loc[candidate] == 2.0:
                        count += 1
                        NGT_df.loc[(NGT_df.index == cell) & (NGT_df['cell_barcode'] == cell_barcode), 'clone_id'] = max_cn
                        NGT_df_homvaf.loc[(NGT_df_homvaf.index == cell) & (NGT_df_homvaf['cell_barcode'] == cell_barcode), 'clone_id'] = max_cn
                        cn_assignment_df.loc[(cn_assignment_df.index == cell) & (cn_assignment_df['cell_barcode'] == cell_barcode), 'clone_id'] = max_cn

    if added_clone == True:
        median_clusters, median_clusters_all = generate_median_cluster_profiles(NGT_df)
        NGT_df, cn_assignment_df, count = merge_large_clusters(NGT_df, NGT_df_homvaf, cn_assignment_df, median_clusters,  count, sample_name)
        NGT_df, cn_assignment_df, count = dissolve_small_clusters(NGT_df, NGT_df_homvaf, cn_assignment_df, median_clusters,  count, sample_name)
    
    print(sample_name, 'number cells changed assignment', count)
    return cn_assignment_df, NGT_df, df_tapestri_clones, added_clone_mut_dict

def main(args):
    analysis_config_yaml = args.c
    with open(analysis_config_yaml, 'r') as f:
        analysis_config = yaml.safe_load(f)

    amplicon_parameters = pd.read_csv(args.w)

    
    dataset = args.d
    local_directory = args.l
    print('Current dataset:', dataset)
    added_clone_mut_dict = defaultdict(str)
    

    cn_assignment_df = pd.read_csv(args.i, index_col=0)
    cn_assignment_df.index.name = 'idx'
    cn_assignment_df = cn_assignment_df.sort_values(by=['idx', 'cell_barcode'])   
    
    df_tapestri_clones = pd.read_csv(args.n, index_col=0)
    
    
    normal_idx = None
    if len(df_tapestri_clones.index[df_tapestri_clones.eq(2.0).all(axis=1)].to_list()) > 0:
        normal_idx =  df_tapestri_clones.index[df_tapestri_clones.eq(2.0).all(axis=1)].to_list()[0]
    
    if normal_idx is not None:
        for c in cn_assignment_df[cn_assignment_df.index == 'RA18_18-11_1']['clone_id'].unique():
            cn_normal_df = cn_assignment_df[cn_assignment_df.index == 'RA18_18-11_1']
            if len(cn_normal_df[cn_normal_df['clone_id'] == c])/cn_normal_df.shape[0] > 0.05:
                cn_assignment_df['clone_id'] = cn_assignment_df['clone_id'].replace(c, normal_idx)
    

    raw_reads = []
    for f in os.listdir(args.r + dataset + '/'):
        raw_read = pd.read_csv(args.r + dataset + '/' + f, sep='\t')
        raw_read.index = [f.split('.')[0]] * len(raw_read)
        raw_reads.append(raw_read)
    raw_reads_df = pd.concat(raw_reads, axis=0)
    raw_reads_df.index.name = 'idx'
    
   


    #HOMDEL SPLIT
    raw_reads_df = raw_reads_df.merge(cn_assignment_df, on=['idx', 'cell_barcode'], how='inner')
    print('raw_reads', raw_reads_df.shape)
    homdel_seps = []
    muts = [c for c in raw_reads_df.columns if c.startswith('AMPL')]
    for c_id in raw_reads_df['clone_id'].unique():
        clone_df = raw_reads_df[raw_reads_df['clone_id'] == c_id]
        mut_idx = 0
        homdels = []
        split_cells = None
        curr_gene = None
        while mut_idx < len(muts):
            curr_mut = muts[mut_idx]
            ampl_gene = amplicon_parameters[amplicon_parameters['amplicon'] == curr_mut]['gene'].values[0]
            if len(clone_df[curr_mut][clone_df[curr_mut] == 0]) > 0.05 * raw_reads_df.shape[0]:
                if split_cells is None:
                    split_cells = clone_df['cell_barcode'][clone_df[curr_mut] == 0]
                    homdels.append(curr_mut)
                    curr_gene = ampl_gene
                else:
                    orig_split_cells = split_cells.copy()
                    split_cells = pd.merge(split_cells, clone_df['cell_barcode'][clone_df[curr_mut] == 0], on=['idx', 'cell_barcode'], how='inner')
                    if len(split_cells) > 0.05 * raw_reads_df.shape[0] and curr_gene == ampl_gene:
                        homdels.append(curr_mut)
                    else:
                        if len(homdels) > 4:
                            homdel_seps.append((c_id, homdels, orig_split_cells, clone_df.shape[0]))
                        homdels = []
                        split_cells = None
                        curr_gene = None
            mut_idx += 1
    
    max_cn = max(list(cn_assignment_df['clone_id'].unique()))
    
    for (x,y,z,a) in homdel_seps:
        for index, row in z.iterrows():
            cn_assignment_df.loc[(cn_assignment_df.index == index) & (cn_assignment_df['cell_barcode'] == row['cell_barcode']), 'clone_id'] = max_cn + 1
        new_clone_profile = df_tapestri_clones.loc[x].copy()
        y = list(set(y).intersection(set([c for c in df_tapestri_clones.columns if c.startswith('AMPL')])))
        new_clone_profile[y] = 0.0
        df_tapestri_clones = pd.concat([df_tapestri_clones, pd.DataFrame([new_clone_profile])])
        df_tapestri_clones.index.values[-1] = max_cn + 1
        print(x, len(y), len(z), a)
        
        max_cn += 1
    
    samples = cn_assignment_df.index.to_series()
    samples = samples.unique()
    
    cravat_df = None
    sample_obj = None
    hdf5_f = None
    for file in os.listdir(local_directory + dataset + '/'):
        if file.startswith(dataset):
            # CRAVAT File
            if file.endswith('.txt'):
                cravat_f = local_directory + dataset + '/' + file
                cravat_df = pd.read_csv(cravat_f, sep='\t', index_col=0, header=[0,1])
            # HDF5 File
            if file.endswith('.h5'):
                hdf5_f = local_directory + dataset + '/' + file
                sample_obj = mio.load(hdf5_f)
                sample_obj_homvaf = mio.load(hdf5_f)
                sample_obj.dna.genotype_variants(min_dp = 8, min_alt_read = 3, assign_low_conf_genotype=False)
                sample_obj_homvaf.dna.genotype_variants(min_dp = 8, het_vaf=5, hom_vaf = 99)# assign_low_conf_genotype=False)
        

    fig, germline_mutations = get_figure_and_germline_mutations(cravat_df, sample_obj, analysis_config, args.s)

    # Initialize NGT dataframe
    NGT_df = sample_obj.dna.get_attribute('NGT', constraint='row')
    NGT_df['idx'] = NGT_df.index.str.split('-').str[1] + '-' + NGT_df.index.str.split('-').str[2]
    NGT_df['cell_barcode'] = NGT_df.index.str.split('-').str[0] 
    NGT_df = NGT_df.reset_index(drop=True)
    NGT_df.set_index('idx', inplace=True)
    NGT_df = pd.merge(NGT_df, cn_assignment_df,on=['idx', 'cell_barcode'], how='inner')
    NGT_df = NGT_df[NGT_df['clone_id'].notna()]



    NGT_df_homvaf = sample_obj_homvaf.dna.get_attribute('NGT', constraint='row')
    NGT_df_homvaf['idx'] = NGT_df_homvaf.index.str.split('-').str[1] + '-' + NGT_df_homvaf.index.str.split('-').str[2]
    NGT_df_homvaf['cell_barcode'] = NGT_df_homvaf.index.str.split('-').str[0] 
    NGT_df_homvaf = NGT_df_homvaf.reset_index(drop=True)
    NGT_df_homvaf.set_index('idx', inplace=True)
    NGT_df_homvaf = pd.merge(NGT_df_homvaf, cn_assignment_df,on=['idx', 'cell_barcode'], how='inner')
    NGT_df_homvaf = NGT_df_homvaf[NGT_df_homvaf['clone_id'].notna()]
 
    
    bin_NGT_df = binarize_NGT_matrix(NGT_df)
    name = ""
    #cn_assignment_df, NGT_df, df_tapestri_clones, added_clone_mut_dict = split_sample_clones(name, bin_NGT_df[germline_mutations + ['cell_barcode', 'clone_id']], NGT_df_homvaf[germline_mutations + ['cell_barcode', 'clone_id']], df_tapestri_clones, cn_assignment_df, added_clone_mut_dict)
    cn_assignment_df, NGT_df, df_tapestri_clones, added_clone_mut_dict = split_sample_clones(name, bin_NGT_df, NGT_df_homvaf, df_tapestri_clones, cn_assignment_df, added_clone_mut_dict)


    #_, median_clusters_all = generate_median_cluster_profiles(bin_NGT_df[germline_mutations + ['cell_barcode', 'clone_id']]) #bin_NGT_df?
    _, median_clusters_all = generate_median_cluster_profiles(bin_NGT_df) #bin_NGT_df?

    #combined_NGT_df, cn_assignment_df, df_tapestri_clones = combine_similar_final_clusters(bin_NGT_df[germline_mutations + ['cell_barcode', 'clone_id']], cn_assignment_df, median_clusters_all, df_tapestri_clones)
    combined_NGT_df, cn_assignment_df, df_tapestri_clones = combine_similar_final_clusters(bin_NGT_df, cn_assignment_df, median_clusters_all, df_tapestri_clones)

    _, median_clusters_all = generate_median_cluster_profiles(combined_NGT_df)
    

    #cn_assignment_df = cn_assignment_df[cn_assignment_df.index != 'RA18_18-11_1']
    n_idx = 0
    new_idx_dict = {}
    for c_id in df_tapestri_clones.index:
        if c_id not in list(cn_assignment_df['clone_id'].unique()):
            df_tapestri_clones = df_tapestri_clones.drop(c_id)
        else:
            new_idx_dict[c_id] = n_idx
            n_idx += 1




    df_tapestri_clones = df_tapestri_clones.reset_index(drop=True)
    cn_assignment_df['clone_id'] = cn_assignment_df['clone_id'].map(new_idx_dict)

    ampl_cols = [c for c in df_tapestri_clones.columns if c.startswith('AMPL')]
    matching_rows = df_tapestri_clones.loc[(df_tapestri_clones[ampl_cols] == len(ampl_cols) * [2.0]).all(axis=1)]

    if len(list(matching_rows.index)) > 0:
        row_index1 = 0
        row_index2 = list(matching_rows.index)[0]
        cn_assignment_df['clone_id'].replace({row_index1: row_index2, row_index2: row_index1}, inplace=True)

        df_tapestri_clones.iloc[row_index1], df_tapestri_clones.iloc[row_index2] = df_tapestri_clones.iloc[row_index2].copy(), df_tapestri_clones.iloc[row_index1].copy()
    
    cn_assignment_df.to_csv(args.a)
    df_tapestri_clones.to_csv(args.p)

    return 0



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', type=str, help='dataset name')
    parser.add_argument('-i', type=str, help='input path to optimal falcon clone assignment output')
    parser.add_argument('-n', type=str, help='input path to optimal falcon clone profile output')
    parser.add_argument('-r', type=str, help='input path to raw read count directory')
    parser.add_argument('-l', type=str, help='input path to hdf5 files directory')
    parser.add_argument('-w', type=str, help='input path to amplicon parameters')
    parser.add_argument('-s', type=str, default=None, help='input path to manually annotated SNVs')
    parser.add_argument('-c', type=str, help='input path to snv selection config parameter file')
    parser.add_argument('-a', type=str, help='output refined clone assignment csv file')
    parser.add_argument('-p', type=str, help='output refined copy number profiles csv file')
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    main(args)
