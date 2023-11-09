#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 5 2022

@author: Palash Sashittal
"""

import gurobipy as gp
import numpy as np
import pandas as pd
import networkx as nx
import itertools
from scipy.stats import betabinom
from collections import defaultdict


# fast version of the constrained dollo solver
class solveFastConstrainedDollo():

    def __init__(self, df_character_matrix, df_total_readcounts = None, df_variant_readcounts = None, force_presence = False,
                 snp_list = [], snv_list = [], annotations = None,
                 k = None, fp = None, fn = None, ado_precision = 15, threads = 1, timelimit = None, verbose = True):
        
        # input character matrix and clustering
        self.df_character_matrix = df_character_matrix
        self.clustering = self.df_character_matrix['cluster_id'].values
        self.A = df_character_matrix.values[:, :-1]
        self.snp_list = snp_list
        self.force_presence = force_presence
        # AKHIL additions
        self.snv_list = snv_list
            
        def mut_replace(x):
            x = x.replace(":", "_").replace("/" , "_").split('_')
            x[2], x[3] = x[3], x[2]
            return "_".join(x)
        
        mapping_dict = list(map(mut_replace, annotations['var_ann'].index))
        
        mapping_dict = dict(zip(annotations['var_ann'].index, mapping_dict))
        # Rename the index using the new_index
        annotations.rename(index=mapping_dict, inplace=True)
        self.mapping_dict = {idx:annotations['var_ann'].loc[idx].values[0] for idx in annotations['var_ann'].index}
    
        self.mutation_list = list(df_character_matrix.columns[:-1])
        
        # input data parameters
        self.ncells = len(self.df_character_matrix)
        self.nclusters = len(self.df_character_matrix['cluster_id'].unique())
        self.nmutations = len(self.df_character_matrix.columns) - 1
        self.k = k
        self.fp = fp
        self.fn = fn        
        
        # read count matrices
        total_reads_mat = df_total_readcounts.values
        alt_reads_mat = df_variant_readcounts.values

        bb_alpha = fp * ado_precision
        bb_beta = (1 - fp) * ado_precision

        presence_coeff_mat = betabinom.logpmf(alt_reads_mat, total_reads_mat, 1, 1)
        absence_coeff_mat = betabinom.logpmf(alt_reads_mat, total_reads_mat, bb_alpha, bb_beta)
        hompresence_coeff_mat = betabinom.logpmf(total_reads_mat - alt_reads_mat, total_reads_mat, bb_alpha, bb_beta)

        plus_state_coeff_mat = np.maximum(presence_coeff_mat, absence_coeff_mat)

        self.cluster_presence_coeff_mat = np.zeros((self.nclusters, self.nmutations))
        self.cluster_absence_coeff_mat = np.zeros((self.nclusters, self.nmutations))
        self.cluster_hompresence_coeff_mat = np.zeros((self.nclusters, self.nmutations))
        self.cluster_plus_state_mat = np.zeros((self.nclusters, self.nmutations))

        for i in range(self.nclusters):
            for j in range(self.nmutations):
                self.cluster_presence_coeff_mat[i,j] = np.sum(presence_coeff_mat[self.clustering == i, j])
                self.cluster_absence_coeff_mat[i,j] = np.sum(absence_coeff_mat[self.clustering == i, j])
                self.cluster_hompresence_coeff_mat[i,j] = np.sum(hompresence_coeff_mat[self.clustering == i, j])        
                self.cluster_plus_state_mat[i,j] = np.sum(plus_state_coeff_mat[self.clustering == i, j])
        
        #self.cluster_plus_state_mat = np.maximum(self.cluster_presence_coeff_mat, self.cluster_absence_coeff_mat)
        
        # gurobi parameters
        self.threads = threads
        self.worklimit = timelimit
        self.verbose = verbose

        # solution
        self.solB = None
        # self.solT_cell = None
        self.solT_mut = None

    def solveSetInclusion(self):

        ncells = self.ncells
        nmutations = self.nmutations
        nclusters = self.nclusters
        
        print(f'n = {ncells}, m = {nmutations}, p = {nclusters}')
        
        clustering = self.clustering
        k = self.k
        snp_list = self.snp_list

        model = gp.Model('solveConstrainedDollo')

        # character matrix variables
        a = model.addVars(nclusters, nmutations, vtype=gp.GRB.BINARY, name='a')
        g = model.addVars(nclusters, nmutations, vtype=gp.GRB.BINARY, name='g')
        c = model.addVars(nclusters, nmutations, k, vtype=gp.GRB.BINARY, name='c')

        # 1+ indicator
        v = model.addVars(nclusters, nmutations, vtype=gp.GRB.CONTINUOUS, lb = 0, ub = 1, name='v')

        # 0 indicator
        w = model.addVars(nclusters, nmutations, vtype=gp.GRB.CONTINUOUS, lb = 0, ub = 1, name='w')

        # column compatibility variables
        y0 = model.addVars(nmutations, nmutations, k, k, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='y0')
        y1 = model.addVars(nmutations, nmutations, k, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='y1')
        y2 = model.addVars(nmutations, nmutations, k, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='y2')
        y3 = model.addVars(nmutations, nmutations, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='y3')

        z0 = model.addVars(nmutations, nmutations, k, k, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='z0')
        z1 = model.addVars(nmutations, nmutations, k, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='z1')
        z2 = model.addVars(nmutations, nmutations, vtype=gp.GRB.CONTINUOUS, lb=0, ub=1, name='z2')

        # mutation can be gained in at most one cluster
        for j in range(nmutations):
            gsum = gp.LinExpr()
            for i in range(nclusters):
                gsum += g[i,j]
            model.addConstr(gsum <= 1)        
        
        # encode one-hot-like constraint on g, a and c
        for i in range(nclusters):
            for j in range(nmutations):
                csum = gp.LinExpr()
                for s in range(k):
                    csum += c[i,j,s]
                model.addConstr(a[i,j] + g[i,j] + csum <= 1)        

        # 1 -set variables
        for i in range(nclusters):
            for j in range(nmutations):
                csum = gp.LinExpr()
                for s in range(k):
                    csum += c[i,j,s]
                # model.addConstr(u[i,j] == g[i,j] + a[i,j] + csum)
                model.addConstr(v[i,j] == a[i,j] + csum)                

        # SNP constraint
        for j in range(nmutations):
            if self.mutation_list[j] in snp_list:
                for i in range(nclusters):
                    model.addConstr(v[i,j] == 1)
                    
        # 0 variable
        for i in range(nclusters):
            for j in range(nmutations):
                csum = gp.LinExpr()
                for s in range(k):
                    csum += c[i,j,s]
                model.addConstr(w[i,j] == 1 - a[i,j] - g[i,j] - csum)        

        # forcing each mutation to be present in the phylogeny
        if self.force_presence:
            for j in range(nmutations):
                msum = gp.LinExpr()
                for i in range(nclusters):
                    csum = gp.LinExpr()
                    for s in range(k):
                        csum += c[i,j,s]
                    msum = a[i,j] + g[i,j] + csum
                model.addConstr(msum >= 1)
                
        # containment constraints
        for j1 in range(nmutations):
            for j2 in range(nmutations):
                if j1 == j2:
                    continue       
                for s1 in range(k):
                    for s2 in range(k):
                        for l in range(nclusters):
                            model.addConstr(y0[j1,j2,s1,s2] <= 1 - c[l,j2,s2] + c[l,j1,s1])

        for j1 in range(nmutations):
            for j2 in range(nmutations):
                if j1 == j2:
                    continue
                for s2 in range(k):
                    for l in range(nclusters):
                        model.addConstr(y1[j1,j2,s2] <= 1 - c[l,j2,s2] + v[l,j1])

        for j1 in range(nmutations):
            for j2 in range(nmutations):
                if j1 == j2:
                    continue                
                for s1 in range(k):
                    for l in range(nclusters):
                        model.addConstr(y2[j1,j2,s1] <= 1 - v[l,j2] + c[l,j1,s1])

        for j1 in range(nmutations):
            for j2 in range(nmutations):
                if j1 == j2:
                    continue
                for i in range(nclusters):
                    model.addConstr(y3[j1,j2] <= 1 - v[i,j2] + v[i,j1])
        
        # disjointness constraints
        for j1 in range(nmutations):
            for j2 in range(nmutations):
                if j1 >= j2:
                    continue                
                for s1 in range(k):
                    for s2 in range(k):
                        for l in range(nclusters):
                            model.addConstr(z0[j1,j2,s1,s2] <= 2 - c[l,j1,s1] - c[l,j2,s2])

                for i in range(nclusters):
                    model.addConstr(z2[j1,j2] <= 2 - v[i,j2] - v[i,j1])

        for j1 in range(nmutations):
            for j2 in range(nmutations):
                if j1 == j2:
                    continue                
                for s2 in range(k):
                    for l in range(nclusters):
                        model.addConstr(z1[j1,j2,s2] <= 2 - c[l,j2,s2] - v[l,j1])        
        
        # perfect phylogeny constraints
        for j1 in range(nmutations):
            for j2 in range(nmutations):
                if j1 >= j2:
                    continue                
                for s1 in range(k):
                    for s2 in range(k):
                        model.addConstr(y0[j1,j2,s1,s2] + y0[j2,j1,s2,s1] + z0[j1,j2,s1,s2] >= 1)
                model.addConstr(y3[j1,j2] + y3[j2,j1] + z2[j1,j2] >= 1)

        for j1 in range(nmutations):
            for j2 in range(nmutations):
                if j1 == j2:
                    continue                
                for s2 in range(k):
                    model.addConstr(y1[j1,j2,s2] + y2[j2,j1,s2] + z1[j1,j2,s2] >= 1)

        # set objective function
        obj_sum = gp.LinExpr()
        for i in range(nclusters):
            for j in range(nmutations):
                obj_sum += self.cluster_plus_state_mat[i,j] * g[i,j]
                obj_sum += self.cluster_presence_coeff_mat[i,j] * a[i,j]
                obj_sum += self.cluster_absence_coeff_mat[i,j] * w[i,j]

                obj_sum += self.cluster_absence_coeff_mat[i,j] * c[i,j,0]
                obj_sum += self.cluster_hompresence_coeff_mat[i,j] * c[i,j,1]        
        
        model.setObjective(obj_sum, gp.GRB.MAXIMIZE)
        
        model.setParam(gp.GRB.Param.Threads, self.threads)
        model.setParam(gp.GRB.Param.Method, 4)
        
        model.setParam(gp.GRB.Param.FeasibilityTol, 1e-6)
        model.setParam(gp.GRB.Param.IntFeasTol, 1e-6)
        model.setParam(gp.GRB.Param.OptimalityTol, 1e-6)
        
        # model.write('test_model_turing.lp')
        
        model.optimize()
        if model.status == gp.GRB.OPTIMAL:

            solg = np.reshape(model.getAttr('x', g).values(), (nclusters, nmutations))
            sola = np.reshape(model.getAttr('x', a).values(), (nclusters, nmutations))
            solc = np.reshape(model.getAttr('x', c).values(), (nclusters, nmutations, k))
            solw = np.reshape(model.getAttr('x', w).values(), (nclusters, nmutations))                        

            df_solb = pd.DataFrame((sola + 2*solc[:,:,0] + 3*solc[:,:,1]).astype(int),
                                   columns=self.df_character_matrix.columns[:-1],
                                   index=range(nclusters), dtype=int)
                        
            self.solB = df_solb
            df_solb_binary = solveFastConstrainedDollo.expand_multi_state_to_binary(df_solb)
            pruned_events = [x for x in df_solb_binary if x not in [f'{y}_1' for y in snp_list]]
            self.solT_mut, self.solT_cell = solveFastConstrainedDollo.generate_perfect_phylogeny(df_solb_binary[pruned_events])
    
    def writeSolution(self, fname):
        if self.solB is not None:
            self.solB.to_csv(fname)
    
    @staticmethod
    def expand_multi_state_to_binary(df_multistate):

        ncells = len(df_multistate)
        binarized_mat = None
        binary_col_dict = {}
        for column in df_multistate.columns:
            max_state = df_multistate[column].max()
            for s in range(1, max_state+1):
                state_col = np.zeros((ncells))
                if s == 1:
                    state_col[df_multistate[column] > 0] = 1
                else:
                    state_col[df_multistate[column] == s] = 1

                binary_col_dict[f'{column}_{s}'] = state_col

        df_binary = pd.DataFrame(binary_col_dict, index = df_multistate.index, dtype=int)
        return df_binary    
    
    @staticmethod
    def generate_perfect_phylogeny(df_binary):

        solT_mut = nx.DiGraph()
        solT_mut.add_node('root')

        solT_cell = nx.DiGraph()
        solT_cell.add_node('root')

        df_binary = df_binary[df_binary.sum().sort_values(ascending=False).index]    

        for cell_id, row in df_binary.iterrows():
            if cell_id == 'root':
                continue

            curr_node = 'root'
            for column in df_binary.columns[row.values == 1]:
                if column in solT_mut[curr_node]:
                    curr_node = column
                else:
                    if column in solT_mut.nodes:
                        raise NameError(f'{column} is being repeated')
                    solT_mut.add_edge(curr_node, column)
                    solT_cell.add_edge(curr_node, column)
                    curr_node = column

            solT_cell.add_edge(curr_node, cell_id)   

        return solT_mut, solT_cell

    def writeDOT(self, dot_file, withcells=True):
        if withcells is True:
            writeTree = self.solT_cell
        else:
            writeTree = self.solT_mut
        
        with open(dot_file, 'w') as output:

            output.write(f'digraph N {{\n')
            output.write(f"\toverlap=\"false\"\n")
            output.write(f"\trankdir=\"TB\"\n")
            
            idx_dict = {}
            idx = 0
            beautify = True
            if writeTree is not None:
                
                edgeset = defaultdict(list)
                for edge in writeTree.edges:
                    edgeset[str(edge[0])].append(edge[1])
                compressed_edgeset = []
                for k,v in edgeset.items():
                    print(k,v)
                    if k == 'root':
                        compressed_edgeset.append(set([k]))
                        print('added root')
                    if len(v) == 1:
                        if str(v[0]).isdigit():
                            compressed_edgeset.append(set([str(v[0])]))
                            print('added terminal state')
                        else:
                            found = False
                            for i, elem in enumerate(compressed_edgeset):
                                if k != 'root':
                                    if k in elem or v[0] in elem:
                                        compressed_edgeset[i] = elem.union(set([k, v[0]]))
                                        print('unioned with', elem)
                                        found = True
                            if found == False:
                                if k != 'root':
                                    compressed_edgeset.append(set([k, v[0]]))
                                else:
                                    compressed_edgeset.append(set([v[0]]))
                    else:
                        for value in v:
                            compressed_edgeset.append(set([str(value)]))
                
                print('edgeset', edgeset)
                print()
                print('compressed_edgeset', compressed_edgeset)
                print()
                if beautify:
                    internal_node_idx = 0
                    node_dict = {}
                    for s in compressed_edgeset:
                        if str(list(s)[0]).isdigit() == False:
                            if str(list(s)[0]) == 'root':
                                output.write(f'\t{idx} [label=\"root\", style=\"bold\"];\n')
                                node_dict['root'] = idx

                            else:
                                output.write(f'\t{idx} [label=\"n{internal_node_idx}\", style=\"bold\"];\n')
                                for i, value in enumerate(s):
                                    node_dict[value] = idx
                                internal_node_idx += 1
                        else:
                            output.write(f'\t{idx} [label=\"')
                            for i, value in enumerate(s):
                                if i == len(s) - 1:
                                    output.write(f'{value}')
                                else:
                                    output.write(f'{value}\\n')
                                node_dict[value] = idx
                            output.write('\", style=\"bold\"];\n')
                        idx_dict[idx] = s
                        idx += 1
                    
                    test = defaultdict(list)
                    for i in self.snp_list:
                        test[self.mapping_dict[i]].append(i)

                    print(node_dict)
                    som_event_dict = {"0": "MISSING", "1": "GAIN", "2": "LOSS", "3": "LOH"}
                    germ_event_dict = {"0": "MISSING", "2": "LOSS", "3": "LOH"}
                    event_color = {"GAIN": "<font color = \"blue\">", "LOSS": "<font color = \"orange\">", "LOH": "<font color = \"black\">", "MISSING": "<font color = \"yellow\">"}
                    for edge in writeTree.edges:
                        if node_dict[str(edge[0])] != node_dict[str(edge[1])]:
                            if str(list(idx_dict[node_dict[str(edge[1])]])[0]).isdigit() == False:
                                output.write(f"\t{node_dict[str(edge[0])]} -> {node_dict[str(edge[1])]} [label=<")
                                for i, value in enumerate(idx_dict[node_dict[str(edge[1])]]):
                                    is_germline = False
                                    color = "<font color = \"black\">"
                                    if value[:-2] in self.snp_list:
                                        color = "<font color = \"green\">"
                                        is_germline = True
                                    elif value[:-2] in self.snv_list:
                                        color = "<font color = \"red\">"
                                    if i == len(idx_dict[node_dict[str(edge[1])]]) - 1:
                                        if is_germline:
                                            output.write(f'{color}{self.mapping_dict[value[:-2]]} </font> {event_color[germ_event_dict[value[-1:]]]} {germ_event_dict[value[-1:]]} </font> ')
                                        else:
                                            output.write(f'{color}{self.mapping_dict[value[:-2]]} </font> {event_color[som_event_dict[value[-1:]]]} {som_event_dict[value[-1:]]}</font>')
                                    else:
                                        if is_germline:
                                            output.write(f'{color}{self.mapping_dict[value[:-2]]} </font> {event_color[germ_event_dict[value[-1:]]]} {germ_event_dict[value[-1:]]} </font><br/>')
                                        else:
                                            output.write(f'{color}{self.mapping_dict[value[:-2]]} </font> {event_color[som_event_dict[value[-1:]]]} {som_event_dict[value[-1:]]}</font><br/>')
                                output.write('>, style=\"bold\"];\n')
                            else: 
                                output.write(f"\t{node_dict[str(edge[0])]} -> {node_dict[str(edge[1])]} [style=\"bold\"];\n")

                else:
                    for node in writeTree.nodes:
                        idx_dict[node] = idx
                        output.write(f'\t{idx} [label=\"{node}\", style=\"bold\"];\n')
                        idx += 1
                    
                    for edge in writeTree.edges:
                        output.write(f"\t{idx_dict[edge[0]]} -> {idx_dict[edge[1]]} [style=\"bold\"];\n")

                output.write(f'}}')
