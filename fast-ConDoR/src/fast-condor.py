#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 5 2021

@author: Palash Sashittal
"""

import pandas as pd
import os
import sys
import argparse
import itertools
import math
import numpy as np
from solveFastConstrainedDollo import solveFastConstrainedDollo
import yaml
# from IPython import embed

def tree_to_newick(T, root=None):
    if root is None:
        roots = list(filter(lambda p: p[1] == 0, T.in_degree()))
        print(roots)
        assert 1 == len(roots)
        root = roots[0][0]
    subgs = []
    while len(T[root]) == 1:
        root = list(T[root])[0]
    for child in T[root]:
        while len(T[child]) == 1:
            child = list(T[child])[0]
        if len(T[child]) > 0:
            child_newick = tree_to_newick(T, root=child)
            if child_newick != '()':
                subgs.append(child_newick)
        else:
            # if child.startswith('s'):
            subgs.append(child)
    # return "(" + ','.join(map(str, subgs)) + ")"
    if len(subgs) == 1:
        return str(subgs[0])
    else:
        return "(" + ','.join(map(str, subgs)) + ")"
    
def main(args):
    df_character_matrix = pd.read_csv(f'{args.i}', index_col = 0)
    df_total_readcounts = pd.read_csv(f'{args.r}', index_col = 0)
    df_variant_readcounts = pd.read_csv(f'{args.v}', index_col = 0)
    
    ncells = len(df_variant_readcounts)
    
    subclonal_mutations = None
    
    if args.subclonal_mutations is not None and os.path.isfile(args.subclonal_mutations):
        with open(args.subclonal_mutations, 'r') as file:
            subclonal_mutations = yaml.safe_load(file)
    else:
        print("[WARNING] manually enforced subclonal mutations file not found. Proceeding with automatic subclonal mutation selection & refinement.")

    cn_profiles = None

    if args.cnp is not None:
        cn_profiles = pd.read_csv(args.cnp, index_col=0)

    snp_list = []
    if args.s is not None:
        with open(args.s, 'r') as inp:
            for line in inp:
                snp_list.append(line.rstrip('\n'))
    
    snv_list = []
    if args.s2 is not None:
        with open(args.s2, 'r') as inp:
            for line in inp:
                snv_list.append(line.rstrip('\n'))
    

    k = args.k
    fp = args.a
    fn = args.b
    ado = args.ado
    
    # filter mutations @HZ: what does this do?
    if args.f:
        df_amplicon = pd.read_csv(args.m, index_col = 0)
        if not "chr" not in df_amplicon.columns:
            if "chrom" in df_amplicon.columns:
                raise IndexError("amplicon metadata file must contain the column named 'chrom' or 'chr'")
            else:
                # rename
                df_amplicon = df_amplicon.rename(columns={"chrom": "chr"})
        
        mutation_list = list(df_variant_readcounts.columns)
        mut_data = []
        for mutation in mutation_list:
            mut_chr = mutation.split('_')[0].lstrip('chr')
            mut_pos = int(mutation.split('_')[1])
            mut_gene = df_amplicon[((df_amplicon['chr'] == mut_chr) &
                                    (df_amplicon['min_pos'] <= mut_pos) &
                                    (df_amplicon['max_pos'] >= mut_pos))]['gene'].values[0]
            mut_data.append([mutation, mut_chr, mut_pos, mut_gene])

        df_mut_meta = pd.DataFrame(mut_data, columns = ['mutation', 'chr', 'pos', 'gene'])

        ## remove chormosome X mutations
        df_mut_meta_selected = df_mut_meta[df_mut_meta['chr'] != 'X']
        df_mut_meta_selected['chr'] = df_mut_meta_selected['chr'].apply(lambda x: int(x))
        df_mut_meta_selected = df_mut_meta_selected.sort_values(['chr', 'pos'])
        sorted_mutation_list = list(df_mut_meta_selected['mutation'])

        df_variant_readcounts = df_variant_readcounts[sorted_mutation_list]
        df_total_readcounts = df_total_readcounts[sorted_mutation_list]
        df_vaf = df_variant_readcounts / df_total_readcounts
        
        presence_theshold = args.pt
        vaf_threshold_1 = args.vt1
        vaf_threshold_0 = args.vt0
        total_read_threshold = args.trt
        missing_fraction_threshold = args.mft    

        df_thresholded = (df_vaf >= vaf_threshold_0).astype(int)
        df_thresholded += (df_vaf >= vaf_threshold_1).astype(int)
        df_thresholded[df_total_readcounts < total_read_threshold] = -1

        selected_mutations = []
        ncells = len(df_thresholded)
        for mutation in df_thresholded:
            hom_presence_fraction = np.sum(df_thresholded[mutation] == 2) / ncells
            absence_fraction = np.sum(df_thresholded[mutation] == 0) / ncells
            het_presence_fraction = np.sum(df_thresholded[mutation] == 1) / ncells
            missing_fraction = np.sum(df_thresholded[mutation] == -1) / ncells

            test_result = True

            if hom_presence_fraction + missing_fraction_threshold > presence_theshold:
                test_result = False

            if het_presence_fraction + missing_fraction_threshold > presence_theshold:
                test_result = False

            if absence_fraction + missing_fraction > presence_theshold:
                test_result = False

            if missing_fraction > missing_fraction_threshold:
                test_result = False

            if test_result:
                selected_mutations.append(mutation)
                # print(mutation, hom_presence_fraction, absence_fraction, het_presence_fraction, missing_fraction, sep='\t')    
    
        df_vaf = df_vaf[selected_mutations]
        df_variant_readcounts = df_variant_readcounts[selected_mutations]
        df_total_readcounts = df_total_readcounts[selected_mutations]
        df_character_matrix = df_character_matrix[selected_mutations + ['cluster_id']]
        snp_list = [x for x in snp_list if x in selected_mutations]
        # AKHIL additions
        snv_list = [x for x in snv_list if x in selected_mutations]
        
    else:
        df_vaf = df_variant_readcounts / df_total_readcounts
    
    # AKHIL additions
        cravat_df = None
        # @HZ 2023-11-14 removed these as they seem to be unnecessary
        # dataset = args.d
        # local_directory = args.c
        # for file in os.listdir(local_directory + dataset + '/'):
        #     if file.startswith(dataset):
        #     # CRAVAT File
        #         if file.endswith('.txt'):
        #             cravat_f = local_directory + dataset + '/' + file
        #             cravat_df = pd.read_csv(cravat_f, sep='\t', index_col=0, header=[0,1])
    # embed() # @HZ 2024-03-25 Gurobi version 11 is not compatible with the current version of the code
    solver = solveFastConstrainedDollo(df_character_matrix, df_total_readcounts=df_total_readcounts,
                                       df_variant_readcounts=df_variant_readcounts,
                                       k=k, fp=fp, fn=fn,
                                       ado_precision = ado, snp_list=snp_list, snv_list=snv_list, annotations = cravat_df, sample=args.d, scr_flag = args.scr, subclonal_mutations=subclonal_mutations, cnp=cn_profiles)

    solver.solveSetInclusion()

    prefix = args.o
    solver.writeSolution(f'{prefix}_B.csv')
    solver.writeDOT(f'{prefix}_tree.dot')
    solver.writeDOT(f'{prefix}_tree_without_cells.dot', withcells=False)
    
    if solver.solT_cell is not None:
        with open(f'{prefix}_tree.newick', 'w') as out:
            out.write(tree_to_newick(solver.solT_cell) + ';')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, help='csv file with mutation matrix and cluster id', required=True)
    parser.add_argument('-r', type=str, help='csv file with total read count matrix', required=True)
    parser.add_argument('-v', type=str, help='csv file with variant read count matrix', required=True)
    parser.add_argument('-s', type=str, help='file containing list of SNPs')
    
    # AKHIL additions
    parser.add_argument('-s2', type=str, help='file containing list of SNVs')
    # parser.add_argument('-c', type=str, help='file path to CRAVAT files', default=None)
    parser.add_argument('-d', type=str, help='current dataset')

    # parser.add_argument('-a', type=float, help='false positive error rate [0.001]', default = 0.001)
    # parser.add_argument('-b', type=float, help='false negative error rate [0.001]', default = 0.001)
    # parser.add_argument('-a', type=float, help='false positive error rate [0.001]', default = 0.025)
    # parser.add_argument('-b', type=float, help='false negative error rate [0.001]', default = 0.025)
    #parser.add_argument('-a', type=float, help='false positive error rate [0.001]', default = 0.1)
    #parser.add_argument('-b', type=float, help='false negative error rate [0.001]', default = 0.1)
    #parser.add_argument('-a', type=float, help='false positive error rate [0.001]', default = 0.025)
    #parser.add_argument('-b', type=float, help='false negative error rate [0.001]', default = 0.025)
    #parser.add_argument('-a', type=float, help='false positive error rate [0.001]', default = 0.01)
    #parser.add_argument('-b', type=float, help='false negative error rate [0.001]', default = 0.01)
    parser.add_argument('-a', type=float, help='false positive error rate [0.001]', default = 0.03)
    parser.add_argument('-b', type=float, help='false negative error rate [0.001]', default = 0.03)

    # parser.add_argument('--ado', type=float, help='precision parameter for ADO', default=15)
    # parser.add_argument('--ado', type=float, help='precision parameter for ADO', default=10)    
    # parser.add_argument('--ado', type=float, help='precision parameter for ADO', default=5)    
    #parser.add_argument('--ado', type=float, help='precision parameter for ADO', default=7.192249441438314)        
    parser.add_argument('--ado', type=float, help='precision parameter for ADO', default=5)        

    parser.add_argument('-k', type=int, help='maximum number of losses for an SNV [2]', default = 2)
    parser.add_argument('-o', type=str, help='output prefix', required=True)
    parser.add_argument('-t', type=int, help='time limit in seconds [1800]', default = 1800)
    
    parser.add_argument('-p', help='force presence of each mutation in the phylogeny? [False]', default=False, action='store_true')
    
    parser.add_argument('-f', help='should the mutations be filtered? [No]', default=False, action='store_true')
    parser.add_argument('-m', type=str, help='metadata file of the amplicons')
    parser.add_argument('--pt', type=float, help='fraction of cells where mutation must be present [0.85]', default=0.85)
    parser.add_argument('--vt0', type=float, help='VAF threshold to call a mutation homozygous [0.75]', default=0.75)
    parser.add_argument('--vt1', type=float, help='VAF thershold to call a mutation absent [0.25]', default=0.25)
    parser.add_argument('--trt', type=int, help='threshold on total number of reads for reliable measurement [10]', default=10)
    parser.add_argument('--mft', type=float, help='fraction of cells where mutation information can be missing [0.2]', default=0.2)
    parser.add_argument('--scr', help='should subclonal refinement be run? [No]', default=False, action='store_true')
    parser.add_argument('--subclonal_mutations', type=str, help='yaml file to dictionary mapping cluster idx to manually selected subclonal SNVs present', default=None)
    parser.add_argument('--cnp', type=str, help='CSV file to copy number profiles of clusters', default=None)

    

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    
    if args.subclonal_mutations is not None and not args.scr:
        raise Exception("--subclonal_mutations requires --scr to be specified.")
    
    if args.i is None and (args.r is None and args.v is None):
        raise Exception("please provide either the binarized mutation matrix, or the total and variant readcount matrices!")
    
    if args.f:
        if args.m is None:
            raise Exception("amplicon metadata file is required for filtering")
    
    main(args)
