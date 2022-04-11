#!/usr/bin/env python3

'''
worm analysis package functions related to reporting tables
'''

import numpy as np
import logging
import matplotlib.pyplot as plt

def find_numTotal(consistent_num_arr, scores_by_conc, scores_by_well, nscores_insystem, num_experiments, num_concentrations, num_replicates):
    if np.all(consistent_num_arr):
        total_nums = np.sum(scores_by_conc[:,:,0,:].reshape(-1, nscores_insystem*num_experiments), axis=1).astype(np.int32)
    else:
        logging.info('Note, since there were not a consistent number of worms in a given well for every day of the experiment, we find the total number of nematodes treated by the following steps'
            + '\nfirst, finding the day with the highest number of worms in each well for each experiment'
            + '\nsecond, summing the different experiments for each well'
            + '\nfinally, summing the triplicate wells for each concentration')
        total_nums = np.sum(np.sum(np.amax(np.sum(scores_by_well, axis=1), axis=1), axis=1).reshape((num_concentrations, num_replicates)), axis=1).astype(np.int32)
    #total_num_in_well = np.amax(np.sum(scores3_by_well, axis=1), axis=1).astype(np.int32) #this is num_wells * num_experiments
    return(total_nums)

def reportTable(df_tablec, df_tabled, rep_exp, expNames, C_day, drug, stage, strain, reportNum, plotIT50, plotLT50, plotIC50, plotLC50, figname_base=''):
    if (reportNum or plotIT50 or plotLT50) and (plotIC50 or plotLC50):
        fig, (ax1, ax2) = plt.subplots(2,1)
        fig.patch.set_visible(False)
        for ax, table_spec in zip([ax1, ax2], [df_tablec, df_tabled]):
            ax.axis("off")
            ax.axis("tight")
            ax.table(cellText = table_spec.values, colLabels=table_spec.columns, rowLabels=table_spec.index, cellLoc = 'center', loc='center')
        filename = 'table_{}_{}_{}_t50repexp_{}_c50day_{}{}.pdf'.format(drug, stage, strain, expNames[rep_exp], C_day, figname_base)

    else:
        fig, ax = plt.subplots()
        fig.patch.set_visible(False)
        ax.axis('off')
        ax.axis('tight')
        if (reportNum or plotIT50 or plotLT50):
            table_spec = df_tablec
            if reportNum and not (plotIT50 or plotLT50):
                filename = 'table_{}_{}_{}{}.pdf'.format(drug, stage, strain, filename_base)
            else:
                filename = 'table_{}_{}_{}_t50repexp_{}{}.pdf'.format(drug, stage, strain, filename_base)
        elif (plotIC50 or plotLC50):
            table_spec = df_tabled
            filename = 'table_{}_{}_{}_c50day_{}{}.pdf'.format(drug, stage, strain, C_day, filename_base)

        ax.table(cellText = table_spec.values, colLabels=table_spec.columns, rowLabels=table_spec.index, cellLoc = 'center', loc='center')

    fig.tight_layout()
    fig.savefig(filename, format='pdf')
    plt.close(fig)
    logging.info("Wrote table(s) of computed values to {}".format(filename))
