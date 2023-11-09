import pandas as pd

PROB_SEQ_ERROR = 0.001
BETABINOM_ALPHA = 0.25
BETABINOM_BETA = 0.25

from scipy.stats import betabinom, binom
import math
# Calculating the optimal sigma

def log_prob_absent(v,t):
    prob = binom.logpmf(v,t,PROB_SEQ_ERROR)
    return prob

def log_prob_present(v,t):
    prob = math.log(betabinom.pmf(v, t, BETABINOM_ALPHA, BETABINOM_BETA))
    return prob

def log_prob_mixed(v,t):
    prob = math.log(0.5 * math.exp(log_prob_absent(v,t)) + (0.5 * math.exp(log_prob_present(v,t))))
    #print(v,t, 'MIXED:', math.exp(prob), math.exp(log_prob_absent(v,t)), math.exp(log_prob_present(v,t)))
    return prob

def compute_LL_solution(BC, result, mutations):
    print(BC)
    print(result)

    def helper(x,v,t):
        if x == 0:
            return log_prob_absent(v,t)
        elif x == 1:
            return log_prob_present(v,t)


    totalLL = 0
    LLs = pd.DataFrame()
    for mutation in mutations:
        Vs = BC['{}_v'.format(mutation)]
        Ts = BC['{}_t'.format(mutation)]
        
        LLs[mutation] = pd.DataFrame(result[mutation]).apply(lambda x: helper(x[mutation], Vs[x.name], Ts[x.name]), axis=1)
        totalLL += sum(LLs[mutation])

    return totalLL


