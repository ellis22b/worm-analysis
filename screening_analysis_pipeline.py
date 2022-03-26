#!/usr/bin/env python3

import numpy as np
import argparse as ap
import logging
import pandas as pd
import datetime
import sys
from matplotlib import rc, rcParams
import worm_analysis_utils
import worm_analysis_scorelines
import worm_analysis_it
import worm_analysis_ic
import worm_analysis_no3
import worm_analysis_pwm
import worm_analysis_stats
import worm_analysis_report

def setup_formatting():
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

def generate_parser():
    parser = ap.ArgumentParser(description = "C. elegans screening analysis and plot maker")
    #arguments about loading data
    parser.add_argument('--toAssess', action='store', dest='toAssess', nargs='+', type=str, required=True, help='the exported files to assess')
    parser.add_argument('--expNames', action='store', dest='expNames', nargs='+', type=str, required=True, help='the names of the experiments passed in toAssess (same order); will be used to annotate plots and tables')
    parser.add_argument('--num_replicates', action='store', dest='num_replicates', type=int, default=3, help='number of replicate wells for each concentration within an experiment. Note this number is assumed to be the same for each concentration')

    #arguments for annotating plots with info about experiments
    parser.add_argument('--drug', action='store', dest='drug', type=str, required=True, help='the drug used in treatment in the assays (ex. ALB, PYR, IVM, or NTZ, etc); will be used to annotate plots and tables; for now this only accepts a single value')
    parser.add_argument('--strain', action='store', dest='strain', type=str, required=True, help='the strain which was treated in the assays (ex. N2 or Hawaii); will be used to annotate plots and tables')
    parser.add_argument('--lifestage', action='store', dest='stage', type=str, required=True, help='the lifestage which was treated in the assays (ex. L1 or L4); will be used to annotate plots and tables')
    parser.add_argument('--concUnits', action='store', dest='concUnits', type=int, required=True, help='use to specify the units of the concentration. 0 is ug/mL, 1 is uM, 2 is M, 3 is mM, 4 is nM, 5 is ng/mL, 6 is mg/mL, 7 is mg/kg, 8 is mg/g')
    parser.add_argument('--molarmass', action='store', dest='molarmass', type=float, required=True, help='the molar mass of the drug used in treatment in the assays. Must be given in the g/mol units')
    parser.add_argument('--densityLiq', action='store', dest='density', type=float, default=1.0, help='the density of the liquid; must be entered in g/mL; default is that of water (1 g/mL)')

    #arguments about what analyses to do
    ## number of worms table
    parser.add_argument('--noReportNum', action='store_false', dest='reportNum', help='add this argument to skip reporting the total number of worms in each concentration (summed across all input experiments), corresponding to Table 1 of the 2017 paper')
    ## line plots
    parser.add_argument('--include_single_exp_plots', action='store_true', dest='isep', help='when plotting the daily motility or lethalty response by concentration, add this argument to plot single experiments as well as the average')
    parser.add_argument('--noLinePlots', action='store_false', dest='linePlots', help='add this argument to skip plotting all line plots for mortality and motility')
    parser.add_argument('--noLinePlotMot', action='store_false', dest='linePlotsMot', help='add this argument to skip plotting motility line plots only')
    parser.add_argument('--noLinePlotMor', action='store_false', dest='linePlotsMor', help='add this argument to skip plotting mortality line plots only')
    ## LC/IC50 plots
    parser.add_argument('--noConc50',  action='store_false', dest='conc50', help='add this argument to skip analysis for and plotting of all inhibitory or lethal concentration analyses (IC50 and LC50)')
    parser.add_argument('--no_plotIC50', action='store_false', dest='plotIC50', help='add this argument to skip plotting the IC50 (3-2-1-0 scoring)')
    parser.add_argument('--no_plotLC50', action='store_false', dest='plotLC50', help='add this argument to skip plotting the LC50 (1-0 scoring)')
    ## LT/IT50 plots
    parser.add_argument('--noTime50', action='store_false', dest='time50', help='add this argument to skip analysis for and plotting of all inhibitory or lethal time analyses (IT50 and LT50)')
    parser.add_argument('--no_plotIT50', action='store_false', dest='plotIT50', help='add this argument to skip plotting the IT50 (3-2-1-0 scoring)')
    parser.add_argument('--no_plotLT50', action='store_false', dest='plotLT50', help='add this argument to skip plotting the LT50 (3-2-1-0 scoring)')
    ## unequal worms prob model

    ## stats

    ## repeating without 3 scores
    parser.add_argument('--no_runNo3', action='store_false', dest='runNo3', help='add this argument to skip running additional analyses where the 3 & 2 scores are combined')

    ##constraining no3 IC parameters
    ###constraining no3 hill slope
    parser.add_argument('--constrainNo3Hill', action='store_true', dest='constrainNo3Hill_bool', help='provide  if an only if you want to constrain the hill slope for motility scoring when re-analyzing without the 3 score (e.g., 2-1-0 scoring) and then use the --no3Hill argument to provide the constrained Hill value')
    parser.add_argument('--no3Hill', action='store', dest='no3Hill', type=float, default=-1.5, required = False, help='the constant/constrained hill value that will be used rather than curve fitting attempting to optimize; recommended use includes when there are fewer than 3 experiments for a given drug/strain/stage treatment combo')
    ###constraining no3 top value

    ###constraining no3 bottom value

    ## spline orders if no fit is possible.
    parser.add_argument('--no3spline_k1', action='store', dest='no3spline_k1', type=int, default=1, required=False, help='the order of the first part (small concentration) of the spline smoothing if there is no IC50 fit for the no3 analysis')
    parser.add_argument('--no3spline_k2', action='store', dest='no3spline_k2', type=int, default=1, required=False, help='the order of the second part (rest of concentrations) of the spline smoothing if there is no IC50 fit for the no3 analysis')

    #arguments specific to LT/IT50 analyses
    parser.add_argument('--representative', action='store', dest='representative', type=int, default=0, help='which number (specifying order/location) of the input files (1, 2, or 3, etc - based on indexing from 1) that is the representative to be used for I/LT50')

    #arguments specific to LC/IC50 analyses
    parser.add_argument('--C_day', action='store', dest='C_day', type=int, default=4, help='the day (index from 1) to assess for inhibitory and/or lethal concentration')
    parser.add_argument('--x0_value', action='store', dest='x0_val', type=float, default=1e-6, help='value to replace the x=0 [] with when transforming x')
    ##constraining parameters
    ###constraining hill slope
    parser.add_argument('--constrainMotHill', action='store_true', dest='constrainMothill_bool', help='provide if and only if you want to constrain the hill slope for motility scoring (e.g., 3-2-1-0 scoring) and then use the --motHill argument to provide the constrained Hill Value')
    parser.add_argument('--motHill', action='store', dest='motHill', type=float, default=-1.5, required = False, help='the constant/constrained hill value that will be used rather than curve fitting attempting to optimize; recommended use includes when there are fewer than 3 experiments for a given drug/strain/stage treatment combo')
    parser.add_argument('--constrainMorHill', action='store_true', dest='constrainMorHill_bool', help='provide if and only if you want to constrain the hill slope for mortality scoring (e.g., 1-0 scoring) and then use the --morHill argument to provide the constrained Hill Value')
    parser.add_argument('--morHill', action='store', dest='morHill', type=float, default=-0.15, required = False, help='the constant/constrained hill value that will be used rather than curve fitting attempting to optimize; recommended use includes when there are fewer than 3 experiments for a given drug/strain/stage treatment combo')
    ### constraining top value

    ### constraining bottom value

    ## spline orders if no fit is possible
    parser.add_argument('--motspline_k1', action='store', dest='motspline_k1', type=int, default=1, required=False, help='the order of the first part (small concentration) of the spline smoothing if there is no IC50 fit for the motility scoring (e.g., 3-2-1-0 scoring) analysis')
    parser.add_argument('--motspline_k2', action='store', dest='motspline_k2', type=int, default=1, required=False, help='the order of the second part (rest of concentrations) of the spline smoothing if there is no IC50 fit for the motility scoring (e.g., 3-2-1-0 scoring) analysis')
    parser.add_argument('--morspline_k1', action='store', dest='morspline_k1', type=int, default=1, required=False, help='the order of the first part (small concentration) of the spline smoothing if there is no LC50 fit for the mortality scoring (e.g., 3-2-1-0 scoring) analysis')
    parser.add_argument('--morspline_k2', action='store', dest='morspline_k2', type=int, default=1, required=False, help='the order of the second part (rest of concentrations) of the spline smoothing if there is no LC50 fit for the mortality scoring (e.g., 3-2-1-0 scoring) analysis')


    #arguments specific to probability model analyses if unequal worm numbers in well
    parser.add_argument('--random_seed', action='store', dest='random_seed', type=int, default=73)

    #arguments specific to stats analyses

    return (parser)

def main():
    setup_formatting()
    parser_args = generate_parser().parse_args()
    analysis_instance = WormAnalysis(parser_args)
    analysis_instance.run_analysis()

class WormAnalysis():
    def __init__(self, parser_args):
        logfile='worm_analysis_{}_{}_{}.txt'.format(drug, stage, strain)
        logging.basicConfig(filename=logfile, level=logging.INFO, filemode='w', format='%(name)s - %(levelname)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
        logging.info("command used to run analysis:\n {}".format(' '.join(sys.argv)))

        self.drug = parser_args.drug
        self.strain = parser_args.strain
        self.stage = parser_args.stage
        self.concUnits = worm_analysis_utils.find_concUnits(parser_args.concUnits)
        self.molarmass = parser_args.molarmass
        self.density = parser_args.density

        self.C_day = parser_args.C_day
        self.x0_value = parser_args.x0_val
        self.representative = parser_args.representative
        self.num_replicates = parser_args.num_replicates
        self.nscore_insystem = 4

        self.reportNum = parser_args.reportNum
        self.linePlots = parser_args.linePlots
        self.linePlotsMot = parser_args.linePlotsMot
        self.linePlotsMor = parser_args.linePlotsMor
        self.include_single_exp_plots = parser_args.isep
        self.time50 = parser_args.time50
        self.time50Mot = parser_args.plotIT50
        self.time50Mor = parser_args.plotLT50
        self.conc50 = parser_args.conc50
        self.conc50Mot = parser_args.plotIC50
        self.conc50Mor = parser_args.plotLC50

        self.load_data(parser_args.toAssess, parser_args.expNames)


    def find_scores_fromdf(self, indf, conc=False, well=False):
        if conc:
            out_arr = np.hstack((np.array(indf.groupby(by = ["Concentration", "Day"]).sum(skipna=True).loc[:, "num0_in_well"]).reshape((self.num_concentrations,self.num_days)),
                                np.array(indf.groupby(by = ["Concentration", "Day"]).sum(skipna=True).loc[:, "num1_in_well"]).reshape((self.num_concentrations,self.num_days)),
                                np.array(indf.groupby(by = ["Concentration", "Day"]).sum(skipna=True).loc[:, "num2_in_well"]).reshape((self.num_concentrations,self.num_days)),
                                np.array(indf.groupby(by = ["Concentration", "Day"]).sum(skipna=True).loc[:, "num3_in_well"]).reshape((self.num_concentrations,self.num_days)))).reshape((self.num_concentrations,self.nscore_insystem,self.num_days))
        elif well:
            out_arr = np.hstack((np.array(indf.groupby(by = ["Well", "Day"]).sum(skipna=True).loc[:, "num0_in_well"]).reshape((self.num_concentrations*self.num_replicates,self.num_days)),
                                np.array(indf.groupby(by = ["Well", "Day"]).sum(skipna=True).loc[:, "num1_in_well"]).reshape((self.num_concentrations*self.num_replicates,self.num_days)),
                                np.array(indf.groupby(by = ["Well", "Day"]).sum(skipna=True).loc[:, "num2_in_well"]).reshape((self.num_concentrations*self.num_replicates,self.num_days)),
                                np.array(indf.groupby(by = ["Well", "Day"]).sum(skipna=True).loc[:, "num3_in_well"]).reshape((self.num_concentrations*self.num_replicates,self.num_days)))).reshape((self.num_concentrations*self.num_replicates,self.nscore_insystem,self.num_days))
        return(out_arr)

    def load_data(self, toAssess, expNames):
        self.num_experiments = len(toAssess)
        self.consistent_num = np.tile(np.nan, self.num_experiments)
        self.expNames = expNames

        for i, file in enumerate(toAssess):
            df = pd.read_csv(file, sep='\t', header=0)
            if i == 0:
                self.uniq_conc = np.unique(df.loc[:,'Concentration'])
                self.mM = worm_analysis_utils.find_mM(self.uniq_conc, self.concUnits, self.molarmass, self.density)
                self.num_concentrations = np.shape(self.uniq_conc)[0]
                self.num_days = np.shape(np.unique(df.loc[:,'Day']))[0]
                assert 0 not in np.unique(np.array(df.loc[:,'Day']).astype(np.int32)), "Remove Day 0 info from input file"

                '''Make a data frame to store computed values in throughout the analysis; index is concentration'''
                units = np.tile(self.concUnits, self.uniq_conc.shape[0])
                index_conctable = [r'\textbf{%s %s}' % (conc,unit) for conc, unit in zip(self.uniq_conc, units)]
                self.df_tablec = pd.DataFrame(index=index_conctable)
                self.df_tablec.index.name='Concentration'

                '''Make a data frame to store computed values in throughout the analysis; index is Day'''
                index_daytable = np.array([ r'\textbf{%s}' % k for k in np.arange(self.num_days+1).astype(str)])
                self.df_tabled = pd.DataFrame(index=index_daytable)
                self.df_tabled.index.name='Day'
                self.df_tabled[r'\textbf{LC50}'] = np.tile('NC', index_daytable.shape[0])
                self.df_tabled[r'\textbf{IC50}'] = np.tile('NC', index_daytable.shape[0])

                '''set some dictionaries for mapping later'''
                self.conc_to_index = {float(c):k for k,c in enumerate(self.uniq_conc)}
                self.index_to_conc = {k:float(c) for k,c in enumerate(self.uniq_conc)}

                self.well_index_to_conc = {w-1: float(c) for w,c in zip(np.array(df.loc[:,'Well']).astype(np.int32)[0:self.num_concentrations*self.num_replicates], np.array(df.loc[:,'Concentration'])[0:self.num_concentrations*self.num_replicates])}
                self.conc_to_well_index = {float(c): w-1 for w,c in zip(np.array(df.loc[:,'Well']).astype(np.int32)[0:self.num_concentrations*self.num_replicates], np.array(df.loc[:,'Concentration'])[0:self.num_concentrations*self.num_replicates])}

                '''start filling arrays with score data'''
                self.scores3_by_conc = self.find_scores_fromdf(df, conc=True)
                self.scores3_by_well = self.find_scores_fromdf(df, well=True)

            else:
                '''continue filling arrays with score data'''
                self.scores3_by_conc = np.stack((self.scores3_by_conc, self.find_scores_fromdf(df, conc=True)), axis=-1)
                self.scores3_by_well = np.stack((self.scores3_by_well, self.find_scores_fromdf(df, well=True)), axis=-1)

            consistent_num = np.nansum(df.loc[:,"numTotal_equal"] == True) == self.num_concentrations*self.num_replicates
            self.consistent_num[i] = consistent_num
        logging.info("Loaded all the data")
        if not np.all(self.consistent_num):
            logging.warning('wells in {} do not have equal numbers of worms in a given well across the length of the experiment'.format(expNames[np.where(~self.consistent_num)]))
            worm_analysis_utils.plot_wormnums(self.scores3_by_well, self.expNames, self.num_days)

    def drive_linePlots(self, figname_base=''):
        if self.linePlots:
            if self.linePlotsMot:
                if self.include_single_exp_plots:
                    for i, exp in enumerate(self.expNames):
                        worm_analysis_scorelines.plotLineCurves(self.motility_index_scores_by_conc[:,;,i], 'isep_motility_{}_{}_{}_{}{}.png'.format(self.drug, self.stage, self.strain, exp, figname_base), r"\textbf{%s %s on %s %s}" %(exp, self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$", "Motility Index Score", 0, self.nscore_insystem-1, 1, np.arange(self.num_days + 1), self.uniq_conc, self.concUnits, self.mM)
                to_plot_mot = worm_analysis_scorelines.get_avg_across_exp(self.motility_index_scores_by_well, self.num_concentrations, self.num_replicates, self.num_experiments, self.num_days)
                worm_analysis_scorelines.plotLineCurves(to_plot_mot,'average_motility_{}_{}_{}{}.png'.format(self.drug, self.stage, self.strain, figname_base), r"\textbf{%s on %s %s}" %(self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$", "Motility Index Score", 0, self.nscore_insystem-1, 1, np.arange(self.num_days + 1), self.uniq_conc, self.concUnits, self.mM)
            if self.linePlotsMor:
                if self.include_single_exp_plots:
                    for i, exp in enumerate(self.expNames):
                        worm_analysis_scorelines.plotLineCurves(self.mortality_scores_by_conc[:,;,i], 'isep_mortality_{}_{}_{}_{}{}.png'.format(self.drug, self.stage, self.strain, exp, figname_base), r"\textbf{%s %s on %s %s}" %(exp, self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$", "\% Alive", 0, 100, 25, np.arange(self.num_days + 1), self.uniq_conc, self.concUnits, self.mM)
                to_plot_mor = worm_analysis_scorelines.get_avg_across_exp(self.mortality_scores_by_well, self.num_concentrations, self.num_replicates, self.num_experiments, self.num_days)
                worm_analysis_scorelines.plotLineCurves(to_plot_mor,'average_mortality_{}_{}_{}{}.png'.format(self.drug, self.stage, self.strain, figname_base), r"\textbf{%s on %s %s}" %(self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$" , "\% Alive", 0, 100, 25, np.arange(self.num_days + 1), self.uniq_conc, self.concUnits, self.mM)

    def drive_survivalTimePlots(self, figname_base=''):
        if self.time50:
            if self.representative  == 0:
                logging.error("Must provide a representative experiment for Survival Analysis using --representative argument")
                exit(1)
            logging.info('Beginning Survival Analysis for {} as the representative experiment'.format(self.expNames[self.representative-1]))
            if self.plotIT50:
                inhibited = np.full((self.num_concentrations, self.num_days+1), 100, dtype=np.float64)
                inhibited, logging_value, IT50 = worm_analysis_it.findSurvivability(inhibited, self.scores3_by_conc, self.representative-1, self.expNames[self.representative-1], self.nscore_insystem, self.consistent_num[self.representative-1], self.num_concentrations, self.num_days, motility=True)
                self.df_tablec[logging_value] = IT50
                logging.info("Added the IT50 values to the table which will be printed later. Column name is 'IT50'"
            if self.plotLT50:
                mortality = np.full((self.num_concentrations, self.num_days+1), 100, dtype=np.float64)
                mortality, logging_value, LT50 = worm_analysis_it.findSurvivability(mortality, self.scores3_by_conc, self.representative-1, self.expNames[self.representative-1], self.nscore_insystem, self.consisten_num[self.representative-1], self.num_concentrations, self.num_days, mortality=True)
                self.df_tablec[logging_value] = LT50
                logging.info("Added the LT50 values to the table which will be printed later. Column name is 'LT50'")
            if self.plotIT50 and self.plotLT50:
                figname_base_spec = '{}_IT50_LT50_' + '{}_{}_{}_{}{}.png'.format(self.drug, self.stage, self.strain, self.expNames[self.representative-1], figname_base)
                self.plotSurvivalTime(inhibited, mortality, figname_base_spec, self.uniq_conc, self.concUnits, self.drug, self.stage, self.strain, np.arange(self.num_days + 1))
            else:
                if self.plotIT50:
                    figname_base_spec = '{}_IT50_' + '{}_{}_{}_{}{}.png'.format(self.drug, self.stage, self.strain, self.expNames[self.representative-1], figname_base)
                    self.plotSurvivalTime(inhibited, np.nan, figname_base_spec, self.uniq_conc, self.concUnits, self.drug, self.stage, self.strain, np.arange(self.num_days + 1), plot_mortality=False)
                if self.plotLT50:
                    logging.warning('Why do you want to plot only the LT50? I suggest plotting IT50 and LT50 together. But here you go')
                    figname_base_spec = '{}_LT50_' + '{}_{}_{}_{}{}.png'.format(self.drug, self.stage, self.strain, self.expNames[self.representative-1], figname_base)
                    self.plotSurvivalTime(np.nan, mortality, figname_base_spec, self.uniq_conc, self.concUnits, self.drug, self.stage, self.strain, np.arange(self.num_days +1), plot_motility=False)
        logging.info('Completed Survival Analysis')

    def drive_inhibitoryConcentrationPlots(self, figname_base=''):

    def run_analysis(self):
        if self.reportNum:
            self.df_tablec[r'\textbf{Number of Nematodes Treated}'] = worm_analysis_report.find_numTotal(self.consistent_num, self.scores3_by_conc, self.scores3_by_well, self.nscore_insystem, self.num_experiments, self.num_concentrations, self.num_replicates)
            logging.info("Added the total number of nematodes treated by concentration to the table which will be printed later. Column name is 'Number of Nematodes Treated'")

        self.motility_index_scores_by_conc, self.motility_index_scores_by_well = worm_analysis_utils.find_scores_driver(self.scores3_by_conc, self.scores3_by_well, self.num_concentrations, self.num_days, self.num_experiments, self.num_replicates, self.nscore_insystem, motility = True)
        self.mortality_scores_by_conc, self.mortality_scores_by_well = worm_analysis_utils.find_scores_driver(self.scores3_by_conc, self.scores3_by_well, self.num_concentrations, self.num_days, self.num_experiments, self.num_replicates, self.nscore_insystem, mortality = True)

        # Drive Line Plots
        self.drive_linePlots()

        # Drive Survival Analysis
        self.drive_survivalTimePlots()



        # Drive Inhibitory concentration analysis

        # Drive no3 score analyses
        ## repeat line plots, survival analysis,inhibitory concentration analysis

        # Drive unequal worms correction

        ## repeat line plots, surivival analysis, inhibitory concentration analysis, no3 analyses also
