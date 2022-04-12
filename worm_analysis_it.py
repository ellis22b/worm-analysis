#!/usr/bin/env python3

'''
worm analysis package IT/LT related functions
'''

import numpy as np
import matplotlib.pyplot as plt
import logging
import worm_analysis_utils

def plotSurvivalTime(inhibited, mortality, figname_base, uniq_conc, concUnits, drug, stage, strain, days_arr, plot_mortality=True, plot_motility=True, ysep=50, ymin=0, ymax=100, ylabel="\% Alive or \% Uninhibited", xlabel='Days', no3 = False, inhibited_og = np.nan):
    for j, conc in enumerate(uniq_conc):
        fig, ax = plt.subplots()

        ax.axhline(50, linestyle=':', color='black')

        title = r"\textbf{%s %s %s on %s %s}" %(conc, concUnits, drug, stage, strain) + r" $\textbf{$\textit{C. elegans}$}$"
        ax = worm_analysis_utils.format_plots(ax, title, xlabel, ylabel, ymin, ymax, ysep, days_arr)

        if no3:
            if plot_motility:
                ax.plot(days_arr, inhibited_og[j], drawstyle = 'steps-post', color = 'orangered', label = r'$\mathrm{Inhibited\;(IT_{50})}$')
                ax.plot(days_arr, inhibited[j], drawstyle='steps-post', color='dodgerblue', label = r'$\mathrm{combined\;3\;\&\;2\;Inhibited\;(IT_{50})}$')
            if plot_mortality:
                ax.plot(days_arr, mortality[j], drawstyle = 'steps-post', color = 'darkviolet', label = r'$\mathrm{Dead-Alive\;(LT_{50})}$')

        else:
            if plot_motility:
                ax.plot(days_arr, inhibited[j], drawstyle = 'steps-post', color = 'orangered', label = r'$\mathrm{Inhibited\;(IT_{50})}$')
            if plot_mortality:
                ax.plot(days_arr, mortality[j], drawstyle = 'steps-post', color = 'darkviolet', label = r'$\mathrm{Dead-Alive\;(LT_{50})}$')


        '''shrink and set legend'''
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.7, box.height*0.65])
        handles, labels = ax.get_legend_handles_labels()
        lgd = ax.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1, 0.94), fontsize=9)

        new_figname = figname_base.format(conc)
        plt.tight_layout()
        fig.savefig(new_figname, bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.close(fig)
        logging.info('Plotted the figure {}'.format(new_figname))


def findSurvivability(toPopulate, toAnalyze, rep_index, expName, nscore_insystem, consistent_num, num_concentrations, num_days, motility=False, mortality=False, no3 = False):
    ''' at each step, only need to decide the number at risk - the number at risk is directly related to the percent survival
        Then correct so that it's monotonically decreasing like graphpad does but in this case with np.minimum.accumulate(a, axis=1)
        Finally, report the day where we cross 50% or U for undefined if we don't cross 50%'''
    #recall toAnalyze = np.zeros((self.num_concentrations, nscore_insystem, self.num_days, self.num_experiments))
    #       toPopulate = np.full((self.num_concentrations, self.num_days+1), 100, dtype=np.float64)
    toAnalyze_expSpec = toAnalyze[:, :, :, rep_index]
    if consistent_num:
        num_total = np.sum(toAnalyze_expSpec[:,:,0], axis=1).reshape((num_concentrations, 1))
    else:
        num_total = np.sum(toAnalyze_expSpec, axis=1)
    if motility:
        num_at_risk = toAnalyze_expSpec[:,nscore_insystem-1,:] #num at risk is only worms of highest score
        num_at_risk_corrected = np.minimum.accumulate(num_at_risk, axis=1)
        if no3:
            logging_value = r'$\mathrm{\textbf{IT50}\;(combined\;3\;\&\;2)}$'
        else:
            logging_value = r'\textbf{IT50}'
    elif mortality:
        num_at_risk = np.sum(toAnalyze_expSpec[:,1:,:], axis=1) #num at risk is any worm of score 1 or greater
        num_at_risk_corrected = np.minimum.accumulate(num_at_risk, axis=1)
        logging_value = r'\textbf{LT50}'
    toPopulate[:, 1:] = num_at_risk_corrected/num_total*100

    '''report day or U'''
    T50 = np.sum(toPopulate <= 50.0, axis=1)
    T50 = (num_days + 1 - T50).astype(str)
    bool_U = T50 == str(num_days + 1)
    T50[bool_U] = 'U'

    return toPopulate, logging_value, T50
