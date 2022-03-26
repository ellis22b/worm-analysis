#!/usr/bin/env python3

'''
worm analysis package IC/LC related functions
'''

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import sem

def transform_X(X, x0_val):
    '''set [0] to x0_val'''
    bool_0 = X == 0
    X[bool_0] += x0_val
    '''return log10 of concentrations'''
    return np.log10(X)

def inhibitorResponse_equation(X, top, bottom, ic50, hillslope=-1):
    '''
    Y=Bottom + (Top-Bottom)/(1+10^((LogIC50-X)*HillSlope))
    where X is log10([]) with x=0 -> x=1e-6 (i.e. call self.transform_X(X))
    '''
    exponent = (np.log10(ic50)- transform_X(X))*hillslope
    equation = bottom + (top-bottom)/(1+(10**exponent))
    return equation

def dose_response_sigmoid(X, X_mid, hill, bottom, top):
    '''
    For plotting found in def 114 here (https://gist.github.com/yannabraham/5f210fed773785d8b638)
    #This is the Hill-Langmuir equation
    #Y = (c+(d-c)/(1+np.exp(b*(np.log(x)-np.log(e)))))
    '''
    Y = (bottom+(top-bottom)/(1+np.exp(hill*(np.log(X)-np.log(X_mid)))))
    return Y

def plotIC(title, figname, concs, averages, sems, uniq_concs, x0_val, concUnits, spline_k1, spline_k2, curve_fit_ic50, curve_fit_hillslope = -1, curve_fit_top= 100, curve_fit_bottom=0, ylabel='\% Uninhibited', ysep=20, ymin=0, ymax=100, fit_possible=True):
    '''Using the Hill-Langmuir equation.'''

    conc_ticks = uniq_concs.copy()
    bool_0 = conc_ticks == 0
    conc_ticks[bool_0] = x0_val
    log_conc_ticks = np.log10(conc_ticks)

    bool_0 = concs == 0
    concs[bool_0] = x0_val
    log_concs = np.log10(concs)

    if fit_possible:
        linspace_x1 = np.log10(np.linspace(x0_val, np.power(10, np.log10(x0_val)+1), 10e2))
        linspace_x2 = np.log10(np.linspace(conc_ticks[1], conc_ticks[-1], 10e4))
        linspace_x1_antilog = np.power(np.tile(10, linspace_x1.shape[0]), linspace_x1)
        linspace_x2_antilog = np.power(np.tile(10, linspace_x2.shape[0]), linspace_x2)

        curve1e = dose_response_sigmoid(linspace_x1_antilog, curve_fit_ic50, curve_fit_hillslope, curve_fit_top, curve_fit_bottom)
        curve2e = dose_response_sigmoid(linspace_x2_antilog, curve_fit_ic50, curve_fit_hillslope, curve_fit_top, curve_fit_bottom)
    else:
        x_s_1 = np.linspace(log_conc_ticks[0], log_conc_ticks[-1], 450)
        spl_1 = make_interp_spline(log_concs, averages, k=spline_k1)
        power_smooth_1 = spl_1(x_s_1)

        x_s_2 = np.linspace(log_concs[1], log_conc_ticks[-1], 300)
        spl_2 = make_interp_spline(log_concs[1:], averages[1:], k=spline_k2)
        power_smooth_2 = spl_2(x_s_2)
        #idx = np.argwhere(np.diff(np.sign(power_smooth_2 - 50))).flatten()
        #print(np.power(10, x_s_2[idx]), flush=True)


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

    #mean values with SEM
    if fit_possible:
        e1 = ax1.errorbar(log_concs[0], averages[0], yerr=sems[0], linestyle='None', marker='o', color='black', capsize=5, clip_on=False)
        e2 = ax2.errorbar(log_concs[1:], averages[1:], yerr=sems[1:], linestyle='None', marker='o', color='black', capsize=5, clip_on=False)
        ax1.plot(linspace_x1, curve1e, c='black', clip_on=False)
        ax2.plot(linspace_x2, curve2e, c='black', clip_on=False)
    else:
        e1 = ax1.errorbar(log_concs[0], averages[0], yerr=sems[0], marker='o', color='black', capsize=5, clip_on=False)
        s1 = ax1.scatter(log_conc_ticks[0], averages[0], c='black')
        s1.set_clip_on(False)
        ax1.plot(x_s_1, power_smooth_1, c='black')

        e2 = ax2.errorbar(log_concs[1:], averages[1:], yerr=sems[1:], marker='o', color='black', capsize=5, clip_on=False)
        s2 = ax2.scatter(log_conc_ticks[1:], averages[1:], c='black')
        ax2.plot(x_s_2, power_smooth_2, c='black', clip_on=False)
        s2.set_clip_on(False)

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
