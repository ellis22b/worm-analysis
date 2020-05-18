#!/usr/bin/env python3

'''
Y=Bottom + (Top-Bottom)/(1+10^((LogIC50-X)*HillSlope))
where X is log10([]) with x=0 -> x=1e-6
'''

import numpy as np
import argparse as ap
from scipy.optimize import curve_fit
from scipy.stats import sem
import logging
import pandas as pd
import datetime
import matplotlib.pyplot as plt
import sys
from matplotlib import rc, rcParams

rc('axes', linewidth=2)
params = {'font.sans-serif': 'Helvetica',
          'font.size': 12,
          'font.weight': 'bold',
          'legend.frameon': False,
          'legend.labelspacing': 1,
          "text.usetex": True,
          'text.latex.preamble': [r'\usepackage{siunitx}',
                           r'\sisetup{detect-all}',
                           r'\usepackage{sansmath}',
                           r'\sansmath']}
rcParams.update(params)

def main():
    parser = generate_parser()
    args = parser.parse_args()
    analysis_instance = WormAnalysis(args.toAssess, args.drug, args.strain, args.stage, args.concUnits, args.reportNum)
    if args.plotLine3 or args.plotLine1:
        analysis_instance.driveLinePlots(args.plotLine3, args.plotLine1, args.isep, args.expNames)
    if args.plotIT50 or args.plotLT50:
        analysis_instance.driveSurvivalTimePlots(args.plotIT50, args.plotLT50, args.rep, args.expNames)
    #if args.plotIC50 or args.plotLC50:
    #    analysis_instance.driveIC(args.plotIC50, args.plotLC50, args.C_day, args.x0_val)


def generate_parser():
    parser = ap.ArgumentParser(description='C elegans analysis')
    parser.add_argument('--toAssess', action='store', dest='toAssess', nargs='+', type=str, required=True, help='the exported files to assess')
    parser.add_argument('--expNames', action='store', dest='expNames', nargs='+', type=str, required=True, help='the names of the experiments passed in toAssess (same order); will be used to annotate plots and tables')
    parser.add_argument('--strain', action='store', dest='strain', type=str, required=True, help='the strain which was treated in the assays (ex. N2 or Hawaii); will be used to annotate plots and tables')
    parser.add_argument('--lifestage', action='store', dest='stage', type=str, required=True, help='the lifestage which was treated in the assays (ex. L1 or L4); will be used to annotate plots and tables')
    parser.add_argument('--drug', action='store', dest='drug', type=str, required=True, help='the drug used in treatment in the assays (ex. ALB, PYR, IVM, or NTZ, etc); will be used to annotate plots and tables; for now this only accepts a single value')
    parser.add_argument('--concUnits', action='store', dest='concUnits', type=int, required=True, help='use to specify the units of the concentration. 0 is ug/mL, 1 is uM, 2 is M, 3 is mM, 4 is nM, 5 is ng/mL, 6 is mg/mL, 7 is mg/kg, 8 is mg/g')
    parser.add_argument('--plotLine3', action='store', dest='plotLine3', type=bool, default=True, help='whether to plot the daily motility response by concentration')
    parser.add_argument('--plotLine1', action='store', dest='plotLine1', type=bool, default=True, help='whether to plot the daily lethality response by concentration')
    parser.add_argument('--include_single_exp_plots', action='store', dest='isep', type=bool, default=False, help='when plotting the daily motility or lethalty response by concentration, whether to plot single experiments as well as the average')
    parser.add_argument('--plotIC50', action='store', dest='plotIC50', type=bool, default=True, help='whether to plot the IC50 (3-2-1-0 scoring)')
    parser.add_argument('--plotLC50', action='store', dest='plotLC50', type=bool, default=True, help='whether to plot the LC50 (1-0 scoring)')
    parser.add_argument('--C_day', action='store', dest='C_day', type=int, default=4, help='the day (index from 1) to assess for inhibitory and/or lethal concentration')
    parser.add_argument('--x0_value', action='store', dest='x0_val', type=float, default=1e-6, help='value to replace the x=0 [] with when transforming x')
    parser.add_argument('--plotIT50', action='store', dest='plotIT50', type=bool, default=True, help='whether to plot the IT50 (3-2-1-0 scoring)')
    parser.add_argument('--plotLT50', action='store', dest='plotLT50', type=bool, default=True, help='whether to plot the LT50 (3-2-1-0 scoring)')
    parser.add_argument('--representative', action='store', dest='rep', type=int, default=0, help='which number (specifying order/location) of the input files (1, 2, or 3, etc - based on indexing from 1) that is the representative to be used for I/LT50')
    parser.add_argument('--reportNum', action='store', dest='reportNum', type=bool, default=True, help='report the total number of worms in each concentration (summed across all input experiments), corresponding to Table 1 of the 2017 paper')

    return parser

class WormAnalysis():
    def __init__(self, toAssess, drug, strain, stage, concUnits, reportNum):
        self.drug = drug
        self.strain = strain
        self.stage = stage
        self.concUnits = concUnits

        self.concUnits_dict = { 0: r"\boldmath$\mu$" + "g/mL",
                                1: r"\boldmath$\mu$" + "M",
                                2: "M",
                                3: "mM",
                                4: "nM",
                                5: "ng/mL",
                                6: "mg/mL",
                                7: "mg/kg",
                                8: "mg/g" }
        self.conc_colors_lo_to_hi = ['black', 'darkorange', 'darkgoldenrod', 'purple', 'limegreen', 'blue']
        self.conc_markers_lo_to_hi = ['s','o', '^', 'v', 'd', 'o']
        self.conc_marker_outline_lo_to_hi = ['black', 'darkorange', 'darkgoldenrod', 'purple', 'limegreen', 'black']
        logfile='worm_analysis_{}_{}_{}.txt'.format(drug, stage, strain)
        logging.basicConfig(filename=logfile, level=logging.INFO, filemode='w', format='%(name)s - %(levelname)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
        logging.info("command used to run analysis:\n {}".format(' '.join(sys.argv)))
        self.load_data(toAssess)
        if reportNum:
            self.find_numTotal()
        self.find_motility_index_score()
        self.find_mortality_score()


    def load_data(self, toAssess):
        self.num_experiments = len(toAssess)
        for i, file in enumerate(toAssess):
            df = pd.read_csv(file, sep='\t', header=0)
            if i == 0:
                self.uniq_conc = np.unique(df.loc[:,'Concentration'])
                self.conc_index = {}
                self.index_to_conc = {}
                for k,c in enumerate(self.uniq_conc):
                    self.conc_index[float(c)] = k
                    self.index_to_conc[k] = float(c)
                self.num_concentrations = len(self.uniq_conc)
                self.num_days = len(np.unique(df.loc[:,'Day']))
                self.scores3_by_well = np.zeros((self.num_concentrations*3, 4, self.num_days, self.num_experiments)) #num_concentrations*3 because concentrations for each experiment should be in triplicate, 4 because of 0-1-2-3 scoring, num_days, and num_experiments
                self.scores3_by_conc = np.zeros((self.num_concentrations, 4, self.num_days, self.num_experiments))
                self.well_index_to_conc = {}
                self.conc_to_well_index = {}
                for well, conc in zip(np.array(df.loc[:, 'Well']).astype(np.int32)[0:self.num_concentrations*3], np.array(df.loc[:,'Concentration'])[0:self.num_concentrations*3]):
                    well_index = well-1
                    self.well_index_to_conc[well_index] = float(conc)
                    self.conc_to_well_index[float(conc)] = well_index
            for j, line in enumerate(open(file)):
                if j != 0:
                    fields=line.strip('\r\n').split('\t')
                    well_index = int(fields[0]) - 1 #index base 0
                    day_index = int(fields[2]) - 1 #index base 0
                    conc_index = self.conc_index[float(fields[1])]
                    #assert conc_index matches self.conc_index[self.well_index_to_conc[well_index]]
                    try:
                        assert conc_index == self.conc_index[self.well_index_to_conc[well_index]]
                    except AssertionError as err:
                        logging.exception("Assertion failed: Inconsistent well to conc. Observed well {} and conc {} when well {} previously conc {}".format(well_index + 1, self.index_to_conc[conc_index], well_index + 1, self.well_index_to_conc[well_index]))
                        raise err

                    self.scores3_by_well[well_index, 0, day_index,i] = int(fields[3]) #store the num in well scoring 0
                    self.scores3_by_conc[conc_index, 0, day_index,i] += int(fields[3])
                    self.scores3_by_well[well_index, 1, day_index,i] = int(fields[4]) #store the num in well scoring 1
                    self.scores3_by_conc[conc_index, 1, day_index,i] += int(fields[4])
                    self.scores3_by_well[well_index, 2, day_index,i] = int(fields[5]) #store the num in well scoring 2
                    self.scores3_by_conc[conc_index, 2, day_index,i] += int(fields[5])
                    self.scores3_by_well[well_index, 3, day_index,i] = int(fields[6]) #store the num in well scoring 3
                    self.scores3_by_conc[conc_index, 3, day_index,i] += int(fields[6])

            num_trues = np.amax(np.unique(df.loc[:,'Well']).astype(np.int32))
            if num_trues != np.sum(df.loc[:,'numTotal_equal']):
                logging.warning('wells in {} do not have equal numbers of worms in a given well across the length of the experiment'.format(file))

    def find_numTotal(self):
        total_nums = np.sum(self.scores3_by_conc[:,:,0,:].reshape(-1, 4*self.num_experiments), axis=1)
        filename = 'table1_equiv_numTreated_{}_{}_{}.txt'.format(self.drug, self.stage, self.strain)
        toWriteTo = open(filename, 'w+')
        toWriteTo.write('Number of nematodes treated\n{} on {} {} _C. elegans_\n'.format(self.drug, self.stage, self.strain))
        for i in range(self.uniq_conc.shape[0]):
            toWriteTo.write('{} {}:\t{}\n'.format(self.uniq_conc[i], self.concUnits_dict[self.concUnits], int(total_nums[i])))
        toWriteTo.close()
        logging.info("Wrote total number of nematodes treated by concentration to the file {}".format(filename))

    def find_motility_index_score(self):
        '''setting motility index score
        Day 0 is automatically a score of 3 for all wells
        weight the number of worms in well of score 0 by 0, the number of worms in well of score 1 by 1, the number of worms in well of score 2 by 2, and the number of worms in well of score 3 by 3
        Sum these weighted numbers and then divide by the total number of worms in well '''
        self.motility_index_scores_by_conc = np.full((self.num_concentrations, self.num_days+1, self.num_experiments), 3.0, dtype=np.float64) #day 0 will be set to 3 automatically, will set days1 onwards at the end
        adjusted_sum = np.sum(self.scores3_by_conc * np.array([0, 1, 2, 3]).reshape((1,-1,1,1)), axis=1) #weight number of worms by the score and sum: 0*num0_in_well + 1*num1_in_well + 2*num2_in_well + 3*num3_in_well
        divided = adjusted_sum/np.sum(self.scores3_by_conc, axis=1) #divide the sum above (adjusted_sum) by the number of worms total in the well (the denominator)
        self.motility_index_scores_by_conc[:, 1:, :] = divided #setting days1 onwards

        self.motility_index_scores_by_well = np.full((self.num_concentrations*3, self.num_days+1, self.num_experiments), 3.0, dtype=np.float64)
        adjusted_sum_by_well = np.sum(self.scores3_by_well * np.array([0, 1, 2, 3]).reshape((1,-1,1,1)), axis=1)
        divided_by_well = adjusted_sum_by_well/np.sum(self.scores3_by_well, axis=1)
        self.motility_index_scores_by_well[:, 1:, :] = divided_by_well

    def find_mortality_score(self):
        '''setting percent alive
        Day 0 is automatically 100 percent alive for all wells
        weight the nubmer of worms in well of score 0 by 0 and weight the number of worms in well of score 1 by 1, the number of worms in well by score 2 by 1, and the nubmer of worms in well of score 3 by 1
        Sum these weighted number to effectively sum the number of worms in the well that are alive
        then divide by the total number of worms in the well and multiply 100 to convert to percentage'''
        self.mortality_scores_by_conc = np.full((self.num_concentrations, self.num_days+1, self.num_experiments), 100.0, dtype=np.float64) #day 0 will be set to 100 automatically, will set days1 onwards at the end
        adjusted_sum = np.sum(self.scores3_by_conc * np.array([0, 1, 1 ,1]).reshape((1, -1, 1, 1)), axis=1) #weight number of worms by the score and sum: 0*num0_in_well + 1*num1_in_well + 1*num2_in_well + 1*num3_in_well where the score 1 denotes alive
        divided_to_percent = adjusted_sum/np.sum(self.scores3_by_conc, axis=1)*100 #divide the sum above (adjusted_sum) by the number of worms total in the well (the denominator), change to percentage
        self.mortality_scores_by_conc[:, 1:, :] = divided_to_percent

        self.mortality_scores_by_well = np.full((self.num_concentrations*3, self.num_days+1, self.num_experiments), 100.0, dtype=np.float64)
        adjusted_sum_by_well = np.sum(self.scores3_by_well * np.array([0, 1, 1, 1]).reshape((1, -1, 1, 1)), axis=1)
        divided_to_percent_by_well = adjusted_sum_by_well/np.sum(self.scores3_by_well, axis=1)*100
        self.mortality_scores_by_well[:, 1:, :] = divided_to_percent_by_well

    def format_plots(self, ax, title, xlabel, ylabel, ymin, ymax, ysep, xticks, format_x = True):

        '''turn off top and right spines'''
        right_side = ax. spines["right"]
        top_side = ax.spines["top"]
        right_side.set_visible(False)
        top_side.set_visible(False)

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


    def plotLineCurves(self, toPlot, figname, title, ylabel, ymin, ymax, ysep, xlabel='Days'):
        '''right now this doesn't include specifying molarity of the highest concentration'''
        fig, ax = plt.subplots()

        days_arr = np.arange(self.num_days + 1)
        ax = self.format_plots(ax, title, xlabel, ylabel, ymin, ymax, ysep, days_arr)

        for i in range(toPlot.shape[0])[::-1]:
            label = "{} {}".format(self.uniq_conc[i], self.concUnits_dict[self.concUnits])
            ax.plot(days_arr, toPlot[i], c=self.conc_colors_lo_to_hi[i], marker=self.conc_markers_lo_to_hi[i], markeredgecolor=self.conc_marker_outline_lo_to_hi[i], label=label, clip_on = False)

        '''shrink and set legend'''
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height*0.75])
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1, 0.5), fontsize=9)

        fig.savefig(figname)
        plt.close(fig)
        logging.info('Plotted the figure {}'.format(figname))

    def driveLinePlots(self, plotLine3, plotLine1, isep, expNames):
        '''Do we want to include the functionality to compute and plot stdev?'''
        if plotLine3:
            if isep:
                for i, exp in enumerate(expNames):
                    self.plotLineCurves(self.motility_index_scores_by_conc[:, :, i], 'isep_motility_{}_{}_{}_{}.png'.format(self.drug, self.stage, self.strain, exp), r"\textbf{%s %s on %s %s}" %(exp, self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$", "Motility Index Score", 0, 3, 1)
            reshaped_mid_by_well = self.motility_index_scores_by_well.reshape((self.num_concentrations, 3, self.num_days+1, self.num_experiments))
            motility_index_across_exp = np.zeros((self.num_concentrations, 3*self.num_experiments, self.num_days+1), dtype=np.float64)
            for j in range(self.num_experiments):
                motility_index_across_exp[:,j*3:(j*3)+3 ,:] = reshaped_mid_by_well[:,:,:,j]
            motility_avg_across_exp = np.mean(motility_index_across_exp, axis=1)
            self.plotLineCurves(motility_avg_across_exp, 'average_motility_{}_{}_{}.png'.format(self.drug, self.stage, self.strain), r"\textbf{%s on %s %s}" %(self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$", "Motility Index Score", 0, 3, 1)

        if plotLine1:
            if isep:
                for i, exp in enumerate(expNames):
                    self.plotLineCurves(self.mortality_scores_by_conc[:, :, i], 'isep_mortality_{}_{}_{}_{}.png'.format(self.drug, self.stage, self.strain, exp), r"\textbf{%s %s on %s %s}" %(exp, self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$" , "\% Alive", 0, 100, 25)
            reshaped_mort_by_well = self.mortality_scores_by_well.reshape((self.num_concentrations, 3, self.num_days+1, self.num_experiments))
            mortality_across_exp = np.zeros((self.num_concentrations, 3*self.num_experiments, self.num_days+1), dtype=np.float64)
            for j in range(self.num_experiments):
                mortality_across_exp[:,j*3:(j*3)+3, :] = reshaped_mort_by_well[:,:,:,j]
            mortality_avg_across_exp = np.mean(mortality_across_exp, axis=1)
            self.plotLineCurves(mortality_avg_across_exp, 'average_mortality_{}_{}_{}.png'.format(self.drug, self.stage, self.strain), r"\textbf{%s on %s %s}" %(self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$" , "\% Alive", 0, 100, 25)

    def plotSurvivalTime(self, inhibited, mortality, figname_base, plot_mortality=True, plot_motility=True, ysep=50, ymin=0, ymax=100, ylabel="\% Alive or \% Uninhibited", xlabel='Days'):
        #LT50 purple, IT50 orange
        for j, conc in enumerate(self.uniq_conc):
            fig, ax = plt.subplots()

            ax.axhline(50, linestyle=':', color='black')

            title = r"\textbf{%s %s %s on %s %s}" %(conc, self.concUnits_dict[self.concUnits], self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$"
            days_arr = np.arange(self.num_days+1)
            ax = self.format_plots(ax, title, xlabel, ylabel, ymin, ymax, ysep, days_arr)

            if plot_motility:
                ax.plot(days_arr, inhibited[j], drawstyle = 'steps-post', color = 'orangered', label = r'$\mathrm{Inhibited\;(IT_{50})}$')
            if plot_mortality:
                ax.plot(days_arr, mortality[j], drawstyle = 'steps-post', color = 'darkviolet', label = r'$\mathrm{Dead-Alive\;(LT_{50})}$')

            '''shrink and set legend'''
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.7, box.height*0.65])
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1, 0.94), fontsize=9)

            new_figname = figname_base.format(conc)
            fig.savefig(new_figname)
            plt.close(fig)
            logging.info('Plotted the figure {}'.format(new_figname))


    def findSurvivability(self, toPopulate, toAnalyze, rep_index, expName, motility=False, mortality=False):
        ''' at each step, only need to decide the number at risk - the number at risk is directly related to the percent survival
            Then correct so that it's monotonically decreasing like graphpad does but in this case with np.minimum.accumulate(a, axis=1)
            Finally, report the day where we cross 50% or U for undefined if we don't cross 50%'''
        #recall toAnalyze = np.zeros((self.num_concentrations, 4, self.num_days, self.num_experiments))
        #       toPopulate = np.full((self.num_concentrations, self.num_days+1), 100, dtype=np.float64)
        toAnalyze_expSpec = toAnalyze[:, :, :, rep_index]
        num_total = np.sum(toAnalyze_expSpec[:,:,0], axis=1).reshape((self.num_concentrations, 1))
        if motility:
            num_at_risk = toAnalyze_expSpec[:,3,:] #num at risk is only worms of score
            num_at_risk_corrected = np.minimum.accumulate(num_at_risk, axis=1)
            logging_value = 'IT50'
        if mortality:
            num_at_risk = np.sum(toAnalyze_expSpec[:,1:,:], axis=1) #num at risk is any worm of score 1, 2, or 3
            num_at_risk_corrected = np.minimum.accumulate(num_at_risk, axis=1)
            logging_value = 'LT50'
        toPopulate[:, 1:] = num_at_risk_corrected/num_total*100

        '''report day or U'''
        filename = 'table2_equiv_{}_{}_{}_{}_{}.txt'.format(logging_value, self.drug, self.stage, self.strain, expName)
        toWriteTo = open(filename, 'w+')
        toWriteTo.write('The {} values for \n{} on {} {} _C. elegans_\n with {} as the representative experiment analyzed\n'.format(logging_value, self.drug, self.stage, self.strain, expName))
        for i in range(self.uniq_conc.shape[0]):
            below50 = np.sum(toPopulate[i] <= 50.0)
            if below50 == 0:
                T50 = 'U'
            else:
                T50 = self.num_days + 1 - below50
            toWriteTo.write('{} {}:\t{}\n'.format(self.uniq_conc[i], self.concUnits_dict[self.concUnits], T50))
        toWriteTo.close()
        logging.info("Wrote the {} values to the file {}".format(logging_value, filename))

        return toPopulate

    def driveSurvivalTimePlots(self, plotIT50, plotLT50, representative, expNames):
        if representative == 0:
            logging.error("Must provide a representative experiment for Survival Analysis using --representative argument")
            exit(1)
        logging.info('Beginning Survival Analysis for {} as the representative experiment'.format(expNames[representative-1]))
        if plotIT50:
            inhibited = np.full((self.num_concentrations, self.num_days+1), 100, dtype=np.float64)
            inhibited = self.findSurvivability(inhibited, self.scores3_by_conc, representative-1, expNames[representative-1], motility=True)
        if plotLT50:
            mortality = np.full((self.num_concentrations, self.num_days+1), 100, dtype=np.float64)
            mortality = self.findSurvivability(mortality, self.scores3_by_conc, representative-1, expNames[representative-1], mortality=True)
        if plotIT50 and plotLT50:
            figname_base = '{}_IT50_LT50_' + '{}_{}_{}_{}.png'.format(self.drug, self.stage, self.strain, expNames[representative-1])
            self.plotSurvivalTime(inhibited, mortality, figname_base)
        else:
            if IT50:
                figname_base = '{}_IT50' + '{}_{}_{}_{}.png'.format(self.drug, self.stage, self.strain, expNames[representative-1])
                self.plotSurvivalTime(inhibited, 0, figname_base, plot_mortality=False)
            else:
                logging.warning('Why do you want to plot only the LT50? I suggest plotting them together. But here you go')
                figname_base = '{}_LT50_' + '{}_{}_{}_{}.png'.format(self.drug, self.stage, self.strain, expNames[representative-1])
                self.plotSurvivalTime(0, mortality, figname_base, plot_motility=False)
        logging.info('Completed Survival Analysis')

    def transform_X(self, X):
        '''set [0] to x0_val'''
        bool_0 = X == 0
        X[bool_0] += self.x0_val
        '''return log10 of concentrations'''
        return np.log10(X)

    def inhibitorResponse_equation(self, X, top, bottom, ic50, hillslope=-1):
        exponent = (np.log10(ic50)- self.transform_X(X))*hillslope
        equation = bottom + (top-bottom)/(1+(10**exponent))
        return equation

    def plotIC(self, title, figname, concs, averages, sems, curve_fit_top, curve_fit_bottom, curve_fit_ic50, curve_fit_hillslope = -1, ylabel='\% Uninhibited', ysep=20, ymin=0, ymax=100):
        linspace_x = np.log10(np.linspace(self.x0_val, np.amax(self.uniq_conc), 1000))
        linspace_x_antilog = np.power(np.tile(10, linspace_x.shape[0]), linspace_x)

        curve = self.inhibitorResponse_equation(linspace_x, curve_fit_top, curve_fit_bottom, curve_fit_ic50, curve_fit_hillslope)

        fig, ax = plt.subplots()

        xlabel='Concentration (%s)' % self.concUnits_dict[self.concUnits]
        ax.axhline(50, linestyle=':', color='black')
        ax = self.format_plots(ax, title, xlabel, ylabel, ymin, ymax, ysep, [], format_x=False)

        #need broken x-axis
        conc_ticks = self.uniq_conc.copy()
        bool_0 = conc_ticks == 0
        conc_ticks[bool_0] += self.x0_val

        #mean values with SEM
        ax.errorbar(concs, averages, yerr=sems, ls='', marker='o', mfc='black', mec='black', clip_on=False)
        ax.plot(linspace_x_antilog, curve, c='black')

        fig.savefig(figname)
        plt.close(fig)

        return 0

    def driveIC(self, plotIC50, plotLC50, C_day, x0_val):
        self.C_day = C_day
        self.x0_val = x0_val





main()
