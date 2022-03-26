#!/usr/bin/env python3

'''
worm analysis package unequal worm number related functions
'''

import numpy as np
import matplotlib as mpl
from matplotlib.text import TextPath
from matplotlib.patches import PathPatch
from matplotlib.font_manager import FontProperties
import worm_analysis_utils

def learn_params(num_concentrations, num_scores, num_days, scores_by_conc):
    '''
    Here we use the number of actually observed worms and the observed scores to learn PCM and PPM parameters for each score for each day for each concentration
    when we generate based on these parameters later, we'll generate the number that should have been observed each day
    '''
    # learn the params for all the concentrations at the same time (params will be 3D instead of 2D)
    PPM_params = np.zeros((num_concentrations, num_scores, num_days + 1))
    PPM_params[:,:,0] = np.tile([0,0,0,1], (num_concentrations, 1)) #on day 0, only scores of 3

    #use self.scores3_by_conc[conc_i, score, day_index, exp_num] to find numerator and denominator
    totals = np.sum(np.sum(scores_by_conc, axis=3), axis=1).reshape((num_concentrations, 1, num_days))
    num_of_score_by_day_and_conc = np.sum(scores_by_conc, axis=3)
    PPM_params[:,:,1:] = num_of_score_by_day_and_conc /totals
    return (PPM_params) #returns (num_conc, num_scores, num_days +1) matrix with entries between 0 and 1

def make_PPM_byWell(PPM_params, conc_to_well_index, num_reps, num_concentrations):
    '''return a matrix that is (num_conc*num_reps, num_scores, num_days+1)'''
    return PPM_byWell


def generate_sequences(PPM_params_byWell, number_to_generate_byWell, num_scores, seed=73):
    np.random.seed(seed)
    sequences = np.full(PPM_params_byWell.shape[0], np.amax(number_to_generate_byWell), np.nan)
    # I think sequences is going to need a day aspect to it too
    to_gen_from = np.arange(0, num_scores+1, 1)
    for


#def PPM_from_PCM(PCM_params, num_seq):
#    PPM_params = PCM_params / num_seq
#    return (PPM_params)

def H_vector_from_PPM(PPM_params):
    H_vector = -1 * np.nansum(PPM_params * np.log2(PPM_params), axis=0)
    return (H_vector)

def R_vector_seqpenality(H_vector, num_seq, addpenalty=False):
    if addpenalty:
        R_vector = 2 - (H_vector + (3 / (2 * np.log(2) * num_seq)))
    else:
        R_vector = 2 - H_vector
    return (R_vector)

def logo_heights(R_vector, PPM_params):
    heights = PPM_params * R_vector
    return (heights)

def logoAt(logo, x, y, yscale=1, ax=None):
    fp = FontProperties(family="Helvetica", weight="bold")
    globscale = 1.35
    logo_alphabet = { "3" : TextPath((-0.305, 0), "3", size=1, prop=fp),
                      "2" : TextPath((-0.384, 0), "2", size=1, prop=fp),
                      "1" : TextPath((-0.35, 0), "1", size=1, prop=fp),
                      "0" : TextPath((-0.366, 0), "0", size=1, prop=fp) }
    logo_colors = {'3': 'orange',
                   '2': 'red',
                   '1': 'blue',
                   '0': 'darkgreen'}

    text = logo_alphabet[logo]
    t = mpl.transforms.Affine2D().scale(1*globscale, yscale*globscale) + mpl.transforms.Affine2D().translate(x,y) + ax.transData
    p = PathPatch(text, lw=0, fc=logo_colors[logo],  transform=t, clip_on = False)
    if ax != None:
        ax.add_artist(p)
    return p

def plotLogo(heights, title):
    fig, ax = plt.subplots()
    ax = worm_analysis_utils.format_plots(ax, title, 'Day', "Information Content (bits)", 0, 2, 1, np.arange(0,8,1))
    #ax.spines["right"].set_visible(False)
    #ax.spines["top"].set_visible(False)
    ax.spines["left"].set_position(("outward", 20))
    ax.spines["bottom"].set_position(("outward", 15))
    #ax.set_ylim(0, 2)
    #ax.set_yticks(np.arange(0, 2.5, 1))
    #ax.set_ylabel("Information Content (bits)")
    #ax.set_xticks(np.arange(0, 8, 1))
    #ax.set_xlabel("Day")
    ax.tick_params(axis="both", pad=15)
    for i in range(np.shape(heights)[1]):
        y=0
        order_heights = np.argsort(heights[:,i])
        heights_sorted = heights[order_heights,i]
        scores = np.arange(np.shape(heights)[0]).astype(str)[order_heights]
        for score, yval in zip(scores, heights_sorted):
            if yval != 0:
                logoAt(score, i, y, yval, ax=ax)
                y += yval
    return(fig)

def run_logoplot(PPM, num_seqs, addpenalty=False):
    '''
    PPM is a concentration specific matrix number of scores X number of days + 1
    '''
    #PPM = PPM_from_PCM(PCM_params, num_seqs)
    H_vector = H_vector_from_PPM(PPM)
    R_vector = R_vector_seqpenality(H_vector, num_seqs, addpenalty=addpenalty)
    heights = logo_heights(R_vector, PPM)
    fig = plotLogo(heights)
    return(fig)

'''
def drive_probModel():
    learn_params()
    for eachconc in concs:
        run_logoplot()
    for eachwell in eachexperimentwells:
        get_concspec_param()
        get_numWorms()
        generate_paths(concspec_params, well_exp_nworms)
    summarizeScores()
    learn_params()

    run_logoplot()
    return_generatedScores()
    return(generated_scores_by_well, generated_scores_by_conc)
'''
