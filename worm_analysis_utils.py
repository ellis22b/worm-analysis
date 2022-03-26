#!/usr/bin/env python3

'''
worm analysis package utility functions
'''

import numpy as np

def find_scores(toAnalyze, ax1, num_days, num_experiments, nscores_insystem, motility = False, mortality = False):
    if motility:
        '''setting motility index score
        Day 0 is automatically a score of 3 for all wells
        weight the number of worms in well of score 0 by 0, the number of worms in well of score 1 by 1, the number of worms in well of score 2 by 2, and the number of worms in well of score 3 by 3
        Sum these weighted numbers and then divide by the total number of worms in well '''
        fill_val = float(nscores_insystem-1)
        sum_multiple = np.arange(nscores_insystem).reshape((1,-1,1,1))
        divided_to_multiple = 1

    if mortality:
        '''setting percent alive
        Day 0 is automatically 100 percent alive for all wells
        weight the nubmer of worms in well of score 0 by 0 and weight the number of worms in well of score 1 by 1, the number of worms in well by score 2 by 1, and the nubmer of worms in well of score 3 by 1
        Sum these weighted number to effectively sum the number of worms in the well that are alive
        then divide by the total number of worms in the well and multiply 100 to convert to percentage'''
        fill_val = 100.0
        sum_multiple = np.hstack((np.array([0]), np.tile(1, nscores_insystem-1))).reshape((1, -1, 1, 1))
        divided_to_multiple = 100

    scores = np.full((ax1, num_days+1, num_experiments), fill_val, dtype=np.float64) #day 0 will be set to 100 automatically, will set days1 onwards at the end
    adjusted_sum = np.sum(toAnalyze * sum_multiple, axis=1) #weight number of worms by the score and sum: 0*num0_in_well + 1*num1_in_well + 1*num2_in_well + 1*num3_in_well where the score 1 denotes alive; or 0*num0_in_well + 1*num1_in_well + 2*num2_in_well + 3*num3_in_well for motility
    divided = adjusted_sum/np.sum(toAnalyze, axis=1)*divided_to_multiple #divide the sum above (adjusted_sum) by the number of worms total in the well (the denominator) and multiply by 100 for percent alive; by 1 for motility
    scores[:, 1:, :] = divided #setting days1 onwards

    return (scores)

def find_scores_driver(scores_by_conc, scores_by_well, num_concentrations, num_days, num_experiments, num_replicates, nscores_insystem, motility = False, mortality = False, conc = True, well = True):
    if conc:
        ax1 = num_concentrations
        toAnalyze = scores_by_conc
        toreturn_scores_by_conc = find_scores(scores_by_conc, ax1, num_days, num_experiments, nscores_insystem, motility = motility, mortality = mortality)
    if well:
        ax1 = num_concentrations * num_replicates
        toreturn_scores_by_well = find_scores(scores_by_well, ax1, num_days, num_experiments, nscores_insystem, motility = motility, mortality = mortality)
    return (toreturn_scores_by_conc, toreturn_scores_by_well)

def find_motility_index_score(scores_by_conc, scores_by_well, num_concentrations, num_days, num_experiments, num_replicates, nscores_insystem):
    '''setting motility index score
    Day 0 is automatically a score of 3 for all wells
    weight the number of worms in well of score 0 by 0, the number of worms in well of score 1 by 1, the number of worms in well of score 2 by 2, and the number of worms in well of score 3 by 3
    Sum these weighted numbers and then divide by the total number of worms in well '''
    motility_index_scores_by_conc = np.full((num_concentrations, num_days+1, num_experiments), float(nscores_insystem-1), dtype=np.float64) #day 0 will be set to 3 automatically, will set days1 onwards at the end
    adjusted_sum = np.sum(scores_by_conc * np.arange(nscores_insystem).reshape((1,-1,1,1)), axis=1) #weight number of worms by the score and sum: 0*num0_in_well + 1*num1_in_well + 2*num2_in_well + 3*num3_in_well
    divided = adjusted_sum/np.sum(scores_by_conc, axis=1) #divide the sum above (adjusted_sum) by the number of worms total in the well (the denominator)
    motility_index_scores_by_conc[:, 1:, :] = divided #setting days1 onwards

    motility_index_scores_by_well = np.full((num_concentrations*num_replicates, num_days+1, num_experiments), float(nscores_insystem-1), dtype=np.float64)
    adjusted_sum_by_well = np.sum(scores_by_well * np.arange(nscores_insystem).reshape((1,-1,1,1)), axis=1)
    divided_by_well = adjusted_sum_by_well/np.sum(scores3_by_well, axis=1)
    motility_index_scores_by_well[:, 1:, :] = divided_by_well

    return (motility_index_scores_by_conc, motility_index_scores_by_well)


def find_mortality_score(scores_by_conc, scores_by_well, num_concentrations, num_days, num_experiments, num_replicates, nscores_insystem):
    '''setting percent alive
    Day 0 is automatically 100 percent alive for all wells
    weight the nubmer of worms in well of score 0 by 0 and weight the number of worms in well of score 1 by 1, the number of worms in well by score 2 by 1, and the nubmer of worms in well of score 3 by 1
    Sum these weighted number to effectively sum the number of worms in the well that are alive
    then divide by the total number of worms in the well and multiply 100 to convert to percentage'''
    mortality_scores_by_conc = np.full((num_concentrations, num_days+1, num_experiments), 100.0, dtype=np.float64) #day 0 will be set to 100 automatically, will set days1 onwards at the end
    adjusted_sum = np.sum(scores_by_conc * np.hstack((np.array([0]), np.tile(1, nscores_insystem-1))).reshape((1, -1, 1, 1)), axis=1) #weight number of worms by the score and sum: 0*num0_in_well + 1*num1_in_well + 1*num2_in_well + 1*num3_in_well where the score 1 denotes alive
    divided_to_percent = adjusted_sum/np.sum(scores_by_conc, axis=1)*100 #divide the sum above (adjusted_sum) by the number of worms total in the well (the denominator), change to percentage
    mortality_scores_by_conc[:, 1:, :] = divided_to_percent

    mortality_scores_by_well = np.full((num_concentrations*num_replicates, num_days+1, num_experiments), 100.0, dtype=np.float64)
    adjusted_sum_by_well = np.sum(scores_by_well * np.hstack((np.array([0]), np.tile(1, nscores_insystem-1))).reshape((1, -1, 1, 1)), axis=1)
    divided_to_percent_by_well = adjusted_sum_by_well/np.sum(scores_by_well, axis=1)*100
    mortality_scores_by_well[:, 1:, :] = divided_to_percent_by_well

    return (mortality_scores_by_conc, mortality_scores_by_well)

def find_concUnits(concUnit_index):
    concUnits_dict = { 0: r"\boldmath$\mu$" + "g/mL",
                            1: r"\boldmath$\mu$" + "M",
                            2: "M",
                            3: "mM",
                            4: "nM",
                            5: "ng/mL",
                            6: "mg/mL",
                            7: "mg/kg",
                            8: "mg/g" }
    return (concUnits_dict[concUnit_index])

def find_mM(uniq_conc, concUnits, molarmass, density):
    '''input: uniq_conc -- vector of unique concentrations
              concUnits -- integer specifing the units of the concentration. 0 is ug/mL, 1 is uM, 2 is M, 3 is mM, 4 is nM, 5 is ng/mL, 6 is mg/mL, 7 is mg/kg, 8 is mg/g
              molarmass --
              density --
       return: mM --
    '''
    '''Questions: 1) what density is reasonable to assume 2) only 1 decimal place??'''
    ''' given x (molarmass) of x g/mole, y max concentration of units corresponding to dictionary key, and z g/mL density
       do the following to find the equivalent mM of the max concentration
       ug/mL -> mM: 1 mole / x g | 1 g / 1e6 ug | y ug / 1 mL | 1e3 mL / 1 L | 1e3 mM / M
       uM -> mM: y uM | 1 mM / 1e3 uM
       M -> mM : y M | 1e3 mM / 1 M
       mM = mM
       nM -> mM: y nM | 1 mM / 1e6 nM
       ng/mL -> mM: 1 mole / x g | 1 g / 1e9 ng | y ng / 1 mL | 1e3 mL / 1L | 1e3 mM / M
       mg/mL -> mM: 1 mole / x g | 1 g / 1e3 mg | y mg / mL | 1e3 mL / 1L | 1e3 mM / M
       mg/kg -> mM: 1 mole / x g | 1 g / 1e3 mg | y mg / kg | 1e3 kg / 1g | z g / 1 mL | 1e3 mL / 1 L | 1e3 mM / M
       mg/g -> mM : 1 mole / x g | 1 g / 1e3 mg | y mg / g | z g / 1 mL | 1e3 mL / 1 L | 1e3 mM / M
       '''
    mM_dict = {find_concUnits(0): 1/molarmass/1e6*np.amax(uniq_conc)*1e3*1e3,
               find_concUnits(1): np.amax(uniq_conc)/1e3,
               find_concUnits(2): np.amax(uniq_conc)*1e3,
               find_concUnits(3): np.amax(uniq_conc),
               find_concUnits(4): np.amax(uniq_conc)/1e6,
               find_concUnits(5): 1/molarmass/1e9*np.amax(uniq_conc)*1e3*1e3,
               find_concUnits(6): 1/molarmass/1e3*np.amax(uniq_conc)*1e3*1e3,
               find_concUnits(7): 1/molarmass/1e3*np.amax(uniq_conc)*1e3*density*1e3*1e3,
               find_concUnits(8): 1/molarmass/1e3*np.amax(uniq_conc)*density*1e3*1e3}

    mM = '{0:.1f}'.format(mM_dict[concUnits])
    return (mM)

def format_plots(ax, title, xlabel, ylabel, ymin, ymax, ysep, xticks, format_x = True):
    '''turn off top and right spines'''
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    '''set title and axis labels'''
    ax.set_title(title, fontsize=12)
    ax.set_xlabel(r'\textbf{%s}' %xlabel, fontsize=10)
    ax.set_ylabel(r'\textbf{%s}' %ylabel, fontsize=10)

    '''set ylim and y ticks/ticklabels'''
    ax.set_ylim(ymin, ymax)

    y_ticklabels = np.arange(ymin, ymax+ysep, ysep)
    ax.set_yticks(y_ticklabels)
    ax.yaxis.set_tick_params(width=2)
    y_ticklabels = y_ticklabels.astype(str)
    for i in range(y_ticklabels.shape[0]):
        y_ticklabels[i] = r'\textbf{%s}'%y_ticklabels[i]
    ax.set_yticklabels(y_ticklabels)

    '''set x ticks/ticklabels'''
    if format_x:
        ax.set_xticks(xticks)
        ax.xaxis.set_tick_params(width=2)
        x_ticklabels = xticks.astype(str)
        for i in range(x_ticklabels.shape[0]):
            x_ticklabels[i] = r'\textbf{%s}'%x_ticklabels[i]
        ax.set_xticklabels(x_ticklabels)

    return ax

def format_concentrations():
    conc_format_dict = {
    "conc_colors_lo_to_hi": ['black', 'darkorange', 'darkgoldenrod', 'purple', 'limegreen', 'blue'],
    "conc_markers_lo_to_hi": ['s','o', '^', 'v', 'd', 'o'],
    "conc_marker_outline_lo_to_hi": ['black', 'darkorange', 'darkgoldenrod', 'purple', 'limegreen', 'black']
    }
    return (conc_format_dict)

def find_totalWorms_byday(scores_mat):
    totals_by_day = np.sum(scores_mat[:,:,:,:], axis=1)
    return (totals_by_day)

def plot_wormnums(scores_by_well, expNames, num_days):
    welltotals_by_day = find_totalWorms_byday(scores_by_well)
    max_totals_by_exp = np.amax(welltotals_by_day, axis=(0, 1))
    for exp_index in np.arange(length(expNames)):
        fig, ax = plt.subplots(figsize=(10,10))
        for well_index in np.arange(welltotals_by_day.shape[0]):
            ax.plot(np.arange(num_days)+1, welltotals_by_day[well_index,:, exp_index], marker='o', c=np.random.rand(3,), label="well {}".format(well_index +1), clip_on=False)
        ax.set_title("Change in worm numbers by well for experiment {}".format(expNames[exp_index]))
        ax.set_ylim(0, max_totals_by_exp[exp_index])
        ax.set_yticks(np.arange(0, max_totals_by_exp[exp_index], 5))
        ax.set_yticklabels(np.arange(0, max_totals_by_exp[exp_index], 5))
        ax.set_ylabel("Number of Worms")
        ax.set_xlabel("Day")
        ax.set_xticks(np.arange(num_days)+1)
        ax.set_xticklabels(np.arange(num_days)+1)
        ax.legend(bbox_to_anchor=[0.99, 0.5], loc='center left', ncol=2)
        plt.tight_layout()
        fig.savefig("well_num_worms_by_day_{}.png".format(expNames[exp_index]))
        plt.close(fig)
