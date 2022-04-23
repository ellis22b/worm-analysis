#!/usr/bin/env python3

'''
worm analysis package stats related functions
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chisquare
import pandas as pd
import logging

def find_num_inhibited(vec_of_scores, nscore_insystem):
    return(np.sum(vec_of_scores[0:nscore_insystem-1]))

def find_opp_ni(vec_of_scores, nscore_insystem):
    return(vec_of_scores[nscore_insystem-1])

def make_inhibited_arr(score_arr, compare_conc_index, compare_day_index, nscore_insystem):
    vec_of_scores = np.sum(score_arr[compare_conc_index, :, compare_day_index, :], axis=1)
    inhibited_arr = np.hstack((find_num_inhibited(vec_of_scores, nscore_insystem), find_opp_ni(vec_of_scores)))
    return (inhibited_arr)

def find_num_mortality(vec_of_scores):
    return(vec_of_scores[0])

def find_opp_nm(vec_of_scores):
    return(vec_of_scores[1:])

def make_mor_arr(score_arr, compare_conc_index, compare_day_index):
    vec_of_scores = np.sum(score_arr[compare_conc_index, :, compare_day_index, :], axis=1)
    mor_arr = np.hstack((find_num_mortality(vec_of_scores), find_opp_nm(vec_of_scores)))
    return (mor_arr)

def chisquare_it(compare_with_vec, expected_vec):
    if np.sum(compare_with_vec >5 ) == 2 and np.sum(expected_vec >5) == 2:
        chisq, pval_raw = chisquare(compare_with_vec, f_exp = expected_vec)
        return(format(pval_raw, ".2e"))
    else:
        return("U")

def compute_num_tests(pval_dfs_dict):
    num_test_counter = 0
    for conc_key in pval_dfs_dict:
        for test_type in ['inh', 'mor']:
            if test_type in pval_dfs_dict[conc_key]:
                pval_df = pval_dfs_dict[conc_key][test_type]
                num_tests_not_performed = (pval_df.values == 'U').sum()
                num_test_counter += (pval_df.size - num_tests_not_performed)
                pval_dfs_dict[conc_key][test_type] = pval_df.replace('U', np.nan)
    return(num_test_counter, pval_dfs_dict)

def multiple_hypothesis_testing_correction(pval_dfs_dict):
    num_tests, pval_dict = compute_num_tests(pval_dfs)
    if num_tests >= 20:
        logging.info('Multiple hypothesis testing correction (Bonferroni) will be performed because there were {} hypothesis tests total'.format(num_tests))
        for conc_key in pval_dict:
            for test_type in ['inh', 'mor']:
                if test_type in pval_dict[conc_key]:
                    pval_df = pval_dict[conc_key][test_type]
                    pval_df_bool = pval_df < (0.05 / num_tests)
                    pval_df_adjpval = pval_df * num_tests
                    pval_dict[conc_key][test_type]['adj.p.val'] = pval_df_adjpval
                    pval_dict[conc_key][test_type]['bool.sig'] = pval_df_bool
    else:
        logging.info('No multiple hypothesis testing correction performed because there were fewer than 20 hypothesis tests')
    return(pval_dict, num_tests)
