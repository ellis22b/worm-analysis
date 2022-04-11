#!/usr/bin/env python3

'''
worm analysis package IC/LC related functions
'''

import numpy as np
from scipy.interpolate import make_interp_spline, BSpline
from scipy.stats import sem
import matplotlib.pyplot as plt
import logging


def find_avg_sem(scores_by_well, C_day, nscores_insystem, num_concentrations, num_replicates, num_experiments, conc_to_index, well_index_to_conc, mortality = False, motility = False):
    totals = np.sum(scores_by_well[:,:,C_day-1,:], axis=1)

    if motility:
        uninhibited = scores_by_well[:,nscores_insystem-1,C_day-1,:]/totals*100
    if mortality:
        uninhibited = np.sum(scores_by_well[:,1:,C_day-1,:], axis=1)/totals*100

    '''find averages and SEMs (note: we plot SEMs because we want to look at the fit of the data to the curve fit model)'''
    toAvg = np.zeros((num_concentrations, num_replicates*num_experiments))
    for i in range(uninhibited.shape[0]):
        toAvg[conc_to_index[well_index_to_conc[i]], ((i%num_replicates) * num_experiments) : (((i%num_replicates)+1) * num_experiments)] = uninhibited[i].copy()

    semdata = sem(toAvg, axis=1)
    avg = np.average(toAvg, axis=1)

    return (uninhibited, avg, semdata)

def set_guess_params(num_exp, observed_averages, notdefhill, notdeftop, notdefbottom, notdefic50, constrainHill = False, motility=False, mortality = False):
    if num_exp < 3 or constrainHill:
        if mortality:
            default_guesses = np.array([np.amax(observed_averages), np.amin(observed_averages), 10**1.5])
        if motility:
            default_guesses = np.array([np.amax(observed_averages), np.amin(observed_averages), 1])

        nondefault_guesses = np.array([notdeftop, notdefbottom, notdefic50])
        nondef_bool = ~np.isnan(nondefault_guesses)
        default_guesses[nondef_bool] = nondefault_guesses[nondef_bool]


        return {'P0': default_guesses,
                'P0_top': default_guesses[0],
                'P0_bottom': default_guesses[1],
                'P0_ic50': default_guesses[2],
                'P0_hill': "constraining"}

    else:
        if mortality:
            default_guesses = np.array([np.amax(observed_averages), np.amin(observed_averages), 10**1.5, -1])
        if motility:
            default_guesses = np.array([np.amax(observed_averages), np.amin(observed_averages), 1, -1])

        nondefault_guesses = np.array([notdeftop, notdefbottom, notdefic50, notdefhill])
        nondef_bool = ~np.isnan(nondefault_guesses)
        default_guesses[nondef_bool] = nondefault_guesses[nondef_bool]

        return {'P0': default_guesses,
                'P0_top': default_guesses[0],
                'P0_bottom': default_guesses[1],
                'P0_hill': default_guesses[2],
                'P0_ic50': default_guesses[3]}


def evaluate_no_fit(top, bottom, ic50, highest_conc):
    no_fit = False
    if top == bottom:
        no_fit = True
        logging.info('No fit possible since top and bottom are equal')
    if bottom > 50:
        no_fit = True
        logging.info('No fit possible since bottom is above 50%')
    if ic50 > highest_conc:
        no_fit = True
        logging.info('No fit possible since reported ic50 is higher than highest concentration assayed')
    return(no_fit)

def find_spline(log_concs, averages, spline_k1, spline_k2):
    x_s_1 = np.linspace(log_concs[0], log_concs[1], 150)
    spl_1 = make_interp_spline(log_concs[:2], averages[:2], k=spline_k1)
    power_smooth_1 = spl_1(x_s_1)

    x_s_2 = np.linspace(log_concs[1], log_concs[-1], 300)
    spl_2 = make_interp_spline(log_concs[1:], averages[1:], k=spline_k2)
    power_smooth_2 = spl_2(x_s_2)
    idx = np.argwhere(np.diff(np.sign(power_smooth_2 - 50))).flatten()
    if len(idx) > 0:
        spline2_ic50 = np.power(10, x_s_2[idx])
    else:
        spline2_ic50 = r'\textgreater' + "{}".format(np.power(10, log_concs[-1]))
    return(x_s_1, power_smooth_1, x_s_2, power_smooth_2, spline2_ic50)

def transform_X(X, x0_val=1e-6):
    '''set [0] to x0_val'''
    bool_0 = X == 0
    X[bool_0] += x0_val
    '''return log10 of concentrations'''
    return np.log10(X)

def dose_response_sigmoid(X, X_mid, hill, bottom, top):
    '''
    For plotting found in def 114 here (https://gist.github.com/yannabraham/5f210fed773785d8b638)
    #This is the Hill-Langmuir equation
    #Y = (c+(d-c)/(1+np.exp(b*(np.log(x)-np.log(e)))))
    '''
    Y = (bottom+(top-bottom)/(1+np.exp(hill*(np.log(X)-np.log(X_mid)))))
    return Y

def plotIC(title, figname, averages, sems, uniq_concs, x0_val, concUnits, spline_k1, spline_k2, popt, constrainedHill = -1, use100_0 = False, useobserved = False, ylabel='\% Uninhibited', ysep=20, ymin=0, ymax=100, fit_possible=True):
    curve_fit_top, curve_fit_bottom, curve_fit_ic50 = popt[0], popt[1], popt[2]
    if use100_0:
        curve_fit_top, curve_fit_bottom = 100, 0
    elif useobserved:
        curve_fit_top, curve_fit_bottom = np.amax(averages), np.amin(averages)
    if len(popt) == 4:
        curve_fit_hillslope = popt[3]
    else:
        curve_fit_hillslope = constrainedHill
    '''Using the Hill-Langmuir equation.'''
    fig = plt.figure(constrained_layout=False)
    widths = [1, 8]
    gs = fig.add_gridspec(1, 2, width_ratios=widths, wspace=0.05)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])

    fig.suptitle(r'\textbf{%s}' %title)
    ax1.set_ylabel(r'\textbf{%s}' %ylabel)
    for ax in [ax1, ax2]:
        ax.axhline(50, linestyle=':', color='black', clip_on=False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.set_ylim(ymin, ymax, ysep)
        ax.set_ylim(ymin, ymax, ysep)
    ax2.spines["left"].set_visible(False)
    ax2.set_yticks([])
    ax2.set_yticklabels([])
    xlabel='Concentration (%s)' % concUnits
    fig.text(0.5, 0.02, r'\textbf{%s}' %xlabel, ha='center', va='center')

    conc_ticks = uniq_concs.copy()
    bool_0 = conc_ticks == 0
    conc_ticks[bool_0] = x0_val
    log_conc_ticks = np.log10(conc_ticks)

    log_concs = log_conc_ticks.copy()

    if fit_possible:
        linspace_x1 = np.log10(np.linspace(x0_val, np.power(10, np.log10(x0_val)+1), 10e2))
        linspace_x2 = np.log10(np.linspace(conc_ticks[1], conc_ticks[-1], 10e4))
        linspace_x1_antilog = np.power(np.tile(10, linspace_x1.shape[0]), linspace_x1)
        linspace_x2_antilog = np.power(np.tile(10, linspace_x2.shape[0]), linspace_x2)

        curve1e = dose_response_sigmoid(linspace_x1_antilog, curve_fit_ic50, curve_fit_hillslope, curve_fit_top, curve_fit_bottom)
        curve2e = dose_response_sigmoid(linspace_x2_antilog, curve_fit_ic50, curve_fit_hillslope, curve_fit_top, curve_fit_bottom)

        e1 = ax1.errorbar(log_concs[0], averages[0], yerr=sems[0], linestyle='None', marker='o', color='black', capsize=5, clip_on=False)
        e2 = ax2.errorbar(log_concs[1:], averages[1:], yerr=sems[1:], linestyle='None', marker='o', color='black', capsize=5, clip_on=False)

        ax1.plot(linspace_x1, curve1e, c='black', clip_on=False)
        ax2.plot(linspace_x2, curve2e, c='black', clip_on=False)

    else:
        x_s_1, power_smooth_1, x_s_2, power_smooth_2, spline2_ic50 = find_spline(log_concs, averages, spline_k1, spline_k2)

        e1 = ax1.errorbar(log_concs[0], averages[0], yerr=sems[0], marker='o', linestyle = 'None', color='black', capsize=5, clip_on=False)
        s1 = ax1.scatter(log_conc_ticks[0], averages[0], marker = 'o', linestyle='None')
        s1.set_clip_on(False)
        ax1.plot(x_s_1, power_smooth_1, c='black')

        e2 = ax2.errorbar(log_concs[1:], averages[1:], yerr=sems[1:], marker='o', linestyle = 'None', color='black', capsize=5, clip_on=False)
        s2 = ax2.scatter(log_conc_ticks[1:], averages[1:], c='black', marker='o', linestyle='None')
        s2.set_clip_on(False)
        ax2.plot(x_s_2, power_smooth_2, c='black', clip_on=False)

    for e in [e1, e2]:
        for b in e[1]:
            b.set_clip_on(False)
    ax1.set_xlim(log_conc_ticks[0], log_conc_ticks[0]+1)
    ax1.set_xticks([log_conc_ticks[0], log_conc_ticks[0]+1])
    lead, power = str(x0_val).split("e-")
    ax1.set_xticklabels([r'$\mathrm{10^{-%s}}$' %str(int(power)), ' '])

    ax2.set_xlim(log_conc_ticks[1], log_conc_ticks[-1])
    ax2.set_xticks(log_conc_ticks[1:])
    ticklabels = np.hstack((conc_ticks[1], conc_ticks[2:].astype(np.int32).astype(str)))
    ax2.set_xticklabels(ticklabels)

    ax1.plot((1,1), (1,1), color='black', clip_on=False) #bottom-left line
    ax2.plot((0,0), (1,1), color='black', clip_on=False) #bottom-right line

    fig.savefig(figname)
    plt.close(fig)
    logging.info('Plotted the figure {}'.format(figname))

    if not fit_possible:
        return(spline2_ic50)
    else:
        return(np.nan)

def set_to_log_value(df_tabled, C_day, bottom, ic50, highest_conc, lowest_nz_conc, motility = False, mortality = False, no3=False):
    if mortality:
        logging_value = r'\textbf{LC50}'
        logging_value2 = 'LC50'
    if motility:
        if no3:
            logging_value = r'\textbf{IC50 (no 3)}'
            logging_value2 = 'IC50 (no 3)'
        else:
            logging_value = r'\textbf{IC50}'
            logging_value2 = 'IC50'
    logging_day = r'\textbf{%s}' % str(C_day)
    if isinstance(ic50, float):
        to_log = '{0:.3f}'.format(ic50)
        if bottom > 50 or ic50 > highest_conc:
            to_log = r'\textgreater' + '{}'.format(highest_conc)
        elif bottom < 50 and ic50 < lowest_nz_conc:
            to_log = r'\textless' + '{}'.format(lowest_nz_conc)
    else:
        to_log = ic50
    df_tabled.loc[logging_day, logging_value] = to_log
    logging.info("Added the %s value to the table which will be printed later. Column name is '%s' and the day/row is '%s'" %(logging_value2, logging_value2, C_day))
    return(df_tabled)
