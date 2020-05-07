#!/usr/bin/env python3

'''
Y=Bottom + (Top-Bottom)/(1+10^((LogIC50-X)*HillSlope))
where X is log10([]) with x=0 -> x=1e-6
'''

import numpy as np
import argparse as ap
from scipy.optimize import curve_fit
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
    analysis_instance = WormAnalysis(args.toAssess, args.logfile, args.drug, args.strain, args.stage, args.concUnits, args.reportNum)
    if args.plotLine3 or args.plotLine1:
        analysis_instance.driveLinePlots(args.plotLine3, args.plotLine1, args.isep, args.expNames)
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
    parser.add_argument('--representative', action='store', dest='rep', type=int, default=0, help='which number (specifying order/location) of the input files (1, 2, or 3, etc) that is the representative to be used for I/LT50')
    parser.add_argument('--reportNum', action='store', dest='reportNum', type=bool, default=True, help='report the total number of worms in each concentration (summed across all input experiments), corresponding to Table 1 of the 2017 paper')
    parser.add_argument('--logfile', action='store', dest='logfile', type=str, default='worm_analysis_{}.txt'.format(datetime.datetime.now()), help='the logfile to store information and any warning messages')

    return parser

class WormAnalysis():
    def __init__(self, toAssess, logfile, drug, strain, stage, concUnits, reportNum):
        self.drug = drug
        self.strain = strain
        self.stage = stage
        self.concUnits = concUnits

        self.concUnits_dict = { 0: r"$\mu$"+"/mL",
                                1: r"$\mu$"+"M",
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
        logging.basicConfig(filename=logfile, level=logging.INFO)
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

    def plotLineCurves(self, toPlot, figname, title, ylabel, ymin, ymax, ysep):
        '''right now this doesn't include specifying molarity of the highest concentration'''
        fig, ax = plt.subplots()
        right_side = ax. spines["right"]
        top_side = ax.spines["top"]
        right_side.set_visible(False)
        top_side.set_visible(False)
        ax.set_title(title)
        days_arr = np.arange(self.num_days + 1)
        for i in range(toPlot.shape[0])[::-1]:
            label = "{} {}".format(self.uniq_conc[i], self.concUnits_dict[self.concUnits])
            ax.plot(days_arr, toPlot[i], c=self.conc_colors_lo_to_hi[i], marker=self.conc_markers_lo_to_hi[i], markeredgecolor=self.conc_marker_outline_lo_to_hi[i], label=label, clip_on = False)
        ax.set_xlabel(r'\textbf{Days}')
        ax.set_ylabel(r'\textbf{%s}' %ylabel)
        ax.set_ylim(ymin, ymax)
        y_ticklabels = np.arange(ymin, ymax+ysep, ysep)
        ax.set_yticks(y_ticklabels)
        y_ticklabels = y_ticklabels.astype(str)
        for i in range(y_ticklabels.shape[0]):
            y_ticklabels[i] = r'\textbf{%s}'%y_ticklabels[i]
        ax.set_yticklabels(y_ticklabels)
        box = ax.get_position()
        #not shrunk version
        #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        #shrunk version
        #ax.set_position([box.x0, box.y0, box.width * 0.8, box.height*0.6])
        #shrunk less version
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height*0.75])
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1, 0.5))
        x_ticklabels = np.arange(self.num_days+1)
        ax.set_xticks(x_ticklabels)
        x_ticklabels = x_ticklabels.astype(str)
        for i in range(x_ticklabels.shape[0]):
            x_ticklabels[i] = r'\textbf{%s}'%x_ticklabels[i]
        ax.set_xticklabels(x_ticklabels)
        ax.xaxis.set_tick_params(width=2)
        ax.yaxis.set_tick_params(width=2)
        fig.savefig(figname)
        plt.close(fig)
        logging.info('Plotted the figure {}'.format(figname))

    def driveLinePlots(self, plotLine3, plotLine1, isep, expNames):
        '''Do we want to include the functionality to compute and plot stdev?'''
        if plotLine3:
            if isep:
                for i, exp in enumerate(expNames):
                    self.plotLineCurves(self.motility_index_scores_by_conc[:, :, i], '{}_isep_motility.png'.format(exp), r"\textbf{%s %s on %s %s}" %(exp, self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$" , "Motility Index Score", 0, 3, 1)
            reshaped_mid_by_well = self.motility_index_scores_by_well.reshape((self.num_concentrations, 3, self.num_days+1, self.num_experiments))
            motility_index_across_exp = np.zeros((self.num_concentrations, 3*self.num_experiments, self.num_days+1), dtype=np.float64)
            for j in range(self.num_experiments):
                motility_index_across_exp[:,j*3:(j*3)+3 ,:] = reshaped_mid_by_well[:,:,:,j]
            motility_avg_across_exp = np.mean(motility_index_across_exp, axis=1)
            self.plotLineCurves(motility_avg_across_exp, 'average_motility_{}_{}_{}.png'.format(self.drug, self.stage, self.strain), r"\textbf{%s on %s %s}" %(self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$" , "Motility Index Score", 0, 3, 1)

        if plotLine1:
            if isep:
                for i, exp in enumerate(expNames):
                    self.plotLineCurves(self.mortality_scores_by_conc[:, :, i], '{}_isep_mortality.png'.format(exp), r"\textbf{%s %s on %s %s}" %(exp, self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$" , "Percent Alive", 0, 100, 25)
            reshaped_mort_by_well = self.mortality_scores_by_well.reshape((self.num_concentrations, 3, self.num_days+1, self.num_experiments))
            mortality_across_exp = np.zeros((self.num_concentrations, 3*self.num_experiments, self.num_days+1), dtype=np.float64)
            for j in range(self.num_experiments):
                mortality_across_exp[:,j*3:(j*3)+3, :] = reshaped_mort_by_well[:,:,:,j]
            mortality_avg_across_exp = np.mean(mortality_across_exp, axis=1)
            self.plotLineCurves(mortality_avg_across_exp, 'average_mortality_{}_{}_{}.png'.format(self.drug, self.stage, self.strain), r"\textbf{%s on %s %s}" %(self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$" , "Percent Alive", 0, 100, 25)


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

    def driveIC(self, plotIC50, plotLC50, C_day, x0_val):
        self.C_day = C_day
        self.x0_val = x0_val





main()
