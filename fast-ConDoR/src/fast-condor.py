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

def tree_to_newick(T, root=None):
    if root is None:
        roots = list(filter(lambda p: p[1] == 0, T.in_degree()))
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
    
    if len(subgs) == 1:
        return str(subgs[0])
    else:
        return "(" + ','.join(map(str, subgs)) + ")"
    
def main(args):
    df_character_matrix = pd.read_csv(f'{args.i}', index_col = 0)
    df_total_readcounts = pd.read_csv(f'{args.r}', index_col = 0)
    df_variant_readcounts = pd.read_csv(f'{args.v}', index_col = 0)
    
    subclonal_mutations = None
    
    if args.subclonal_mutations is not None and os.path.isfile(args.subclonal_mutations):
        with open(args.subclonal_mutations, 'r') as file:
            subclonal_mutations = yaml.safe_load(file)
    else:
        print("[WARNING] subclonal mutations file not found. Proceeding without automatic subclonal mutation selection & refinement.")

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
    
    
    solver = solveFastConstrainedDollo(df_character_matrix, df_total_readcounts=df_total_readcounts,
                                       df_variant_readcounts=df_variant_readcounts,
                                       k=k, fp=fp, fn=fn,
                                       ado_precision = ado, snp_list=snp_list, snv_list=snv_list, sample=args.d, scr_flag = args.scr, subclonal_mutations=subclonal_mutations, cnp=cn_profiles)

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
    parser.add_argument('-s2', type=str, help='file containing list of SNVs')
    parser.add_argument('-c', type=str, help='file path to CRAVAT files')
    parser.add_argument('-d', type=str, help='current dataset')
    parser.add_argument('-a', type=float, help='false positive error rate [0.001]', default = 0.03)
    parser.add_argument('-b', type=float, help='false negative error rate [0.001]', default = 0.03)
    parser.add_argument('--ado', type=float, help='precision parameter for ADO', default=5)        
    parser.add_argument('-k', type=int, help='maximum number of losses for an SNV [2]', default = 2)
    parser.add_argument('-o', type=str, help='output prefix', required=True)
    parser.add_argument('-t', type=int, help='time limit in seconds [1800]', default = 1800)
    parser.add_argument('-p', help='force presence of each mutation in the phylogeny? [False]', default=False, action='store_true')
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
    
    main(args)
