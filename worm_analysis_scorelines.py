#!/usr/bin/env python3

import worm_analysis_utils

'''
worm analysis package score line plot functions
'''

def get_avg_across_exp(to_average, num_concentrations, num_replicates, num_experiments, num_days):
    reshaped_by_well = to_average.reshape((num_concentrations, num_replicates, num_days+1, num_experiments))
    across_exp = np.zeros((num_concentrations, num_replicates*num_experiments, num_days+1), dtype=np.float64)
    for j in range(num_experiments):
        across_exp[:,j*num_replicates:(j*num_replicates)+num_replicates,:] = reshaped_by_well[:,:,:,j]
    avg_across_exp = np.mean(across_exp, axis=1)
    return(avg_across_exp)

def plotLineCurves(toPlot, figname, title, ylabel, ymin, ymax, ysep, days_arr, uniq_concs, concUnits, mM, xlabel='Days'):
    '''Specifies molarity of the highest concentration'''
    fig, ax = plt.subplots()

    ax = worm_analysis_utils.format_plots(ax, title, xlabel, ylabel, ymin, ymax, ysep, days_arr)

    for i in range(toPlot.shape[0])[::-1]:
        if i == toPlot.shape[0]-1 and (concUnits in [worm_analysis_utils.find_concUnits(i) for i in [0,5,6,7,8]]):
            label = "{} {}\n= {} mM".format(uniq_concs[i], concUnits, mM)
        else:
            label = "{} {}".format(uniq_concs[i], concUnits)
        ax.plot(days_arr, toPlot[i], c=worm_analysis_utils.format_concentrations()["conc_colors_lo_to_hi"][i], marker=worm_analysis_utils.format_concentrations()["conc_markers_lo_to_hi"][i], markeredgecolor=worm_analysis_utils.format_concentrations()["conc_marker_outline_lo_to_hi"][i], label=label, clip_on = False)

    '''shrink and set legend'''
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height*0.75])
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1, 0.5), fontsize=9)

    fig.savefig(figname)
    plt.close(fig)
    logging.info('Plotted the figure {}'.format(figname))
