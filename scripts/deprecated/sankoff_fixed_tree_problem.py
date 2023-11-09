import numpy as np
import pandas as pd


'''
Givens: 
    n - # of clones
    k - # of cn_states
    cn_profiles - dataframe of leaf profiles
    T - fixed tree phylogeny
Return:
    Labels(T) - copy number profile labels of internal nodes
'''


def small_parsimony(cn_profiles, T):

    DP = np.zeroes((T.size, k, cn_profiles.shape)) # T.size represents the number of nodes, cn_profiles.shape represents the number of amplicons

    for each leaf node in T:
        for each cn_state:
            for each amplicon:
                DP[leaf node][cn_state][amplicon] = 0 if cn_profiles[leaf node][amplicon] == cn_state else infinity

    # Traverse the tree in post-order (bottom-up)
        for each internal node in T (in post-order):
            for each cn_state:
                for each amplicon
                    min_cost = infinity
                    
                    for each possible_state:
                        cost = DP[left child node][possible_state][amplicon] + DP[right child node][cn_state][amplicon] + dist(possible_state, cn_state)
                        if cost < min_cost:
                            min_cost = cost
                    
                    DP[internal node][cn_state][amplicon] = min_cost

    min_parsimony_score = min(DP[T.root]) #find the minimum value for each amplicon

    # Backtrack to find the copy number states at internal nodes
    Backtrack(T.root, cn_profiles)

    '''
    For recursive backtracking, find the minimum cost path backward using the DP table to find the labellings for the cn states at all internal nodes
    '''


