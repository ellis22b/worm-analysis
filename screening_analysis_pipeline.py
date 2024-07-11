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
from scipy.optimize import curve_fit

def setup_formatting():
    rc('axes', linewidth=2)
    params = {"text.usetex": True,
              "font.sans-serif": 'Helvetica',
              "font.size": 12,
              "font.weight": 'bold',
              "legend.frameon": False,
              "legend.labelspacing": 1,
              "text.latex.preamble": [r'\usepackage{siunitx}',
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
    parser.add_argument('--includeSEM', action='store_true', dest='isem', help='when plotting the daily motility or lethality response by concentration, add this argument to plot SEM bars')
    parser.add_argument('--includeSTDEV', action='store_true', dest='istdev', help='when plotting the daily motility or lethality response by concentration, add this argument to plot STDEV bars')
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
    parser.add_argument('--no_runcombine3_2', action='store_false', dest='runNo3', help='add this argument to skip running additional analyses where the 3 & 2 scores are combined')

    ##constraining no3 IC parameters
    ###constraining no3 hill slope
    parser.add_argument('--constrain_combine32Hill', action='store_true', dest='constrainNo3Hill_bool', help='provide  if an only if you want to constrain the hill slope for motility scoring when re-analyzing without the 3 score (e.g., 2-1-0 scoring) and then use the --no3Hill argument to provide the constrained Hill value')
    parser.add_argument('--combine32Hill', action='store', dest='no3Hill', type=float, default=np.nan, required = False, help='the constant/constrained hill value that will be used rather than curve fitting attempting to optimize; required use when there are fewer than 3 experiments for a given drug/strain/stage treatment combo')
    parser.add_argument('--combine32notDefHill', action='store', dest='no3notDefHill', type=float, default=np.nan, required = False, help='not default Hill slope initial guess parameter for IC50 analysis with no3 scoring.')
    ###constraining no3 top or bottom value
    parser.add_argument('--combine32notDefTop', action='store', dest='no3notDefTop', type=float, default=np.nan, required=False, help='not default Top initial guess parameter for IC50 analysis with no3 scoring. Default is the observed value in the dataset')
    parser.add_argument('--combine32notDefBottom', action='store', dest='no3notDefBottom', type=float, default=np.nan, required=False, help='not default Bottom initial guess parameter for IC50 analysis with no3 scoring. Default is the observed value in the dataset')
    parser.add_argument('--plotcombine32TopBottom', action='store_true', dest='plotno3TopBottom', help='use this flag if you want to plot with the true/observed top value rather than the fit top for no3 mortality scoring')
    parser.add_argument('--plotcombine32HundredZero', action='store_true', dest='plotno3HundredZero', help='use this flag if you want to plot with a top of 100 and a bottom of 0 no matter the true/observed top/bottom values or the curve fit values for no3 mortality scoring; do not use this flag with --plotno3TopBottom')

    ###guess no3 ic50
    parser.add_argument('--combine32notDefic50', action='store', dest='no3notDefic50', type=float, default=np.nan, required=False, help='not default ic50 initial guess parameter for IC50 analysis with no3 scoring.')

    ## spline orders if no fit is possible.
    parser.add_argument('--combine32spline_k1', action='store', dest='no3spline_k1', type=int, default=1, required=False, help='the order of the first part (small concentration) of the spline smoothing if there is no IC50 fit for the no3 analysis')
    parser.add_argument('--combine32spline_k2', action='store', dest='no3spline_k2', type=int, default=1, required=False, help='the order of the second part (rest of concentrations) of the spline smoothing if there is no IC50 fit for the no3 analysis')

    #arguments specific to LT/IT50 analyses
    parser.add_argument('--representative', action='store', dest='representative', type=int, default=0, help='which number (specifying order/location) of the input files (1, 2, or 3, etc - based on indexing from 1) that is the representative to be used for I/LT50')

    #arguments specific to LC/IC50 analyses
    parser.add_argument('--C_day', action='store', dest='C_day', type=int, default=4, help='the day (index from 1) to assess for inhibitory and/or lethal concentration')
    parser.add_argument('--x0_value', action='store', dest='x0_val', type=float, default=1e-6, help='value to replace the x=0 [] with when transforming x')
    parser.add_argument('--motReturnliteral50', action='store_true', dest='motreturnliteral50', help='use this flag to correct the motility ic50 reporting if your ic50 fit reports the 50% decrease from top to bottom but not where the line crosses 50%')
    parser.add_argument('--morReturnliteral50', action='store_true', dest='morreturnliteral50', help='use this flag to correct the mortality lc50 reporting if your ic50 fit reports the 50% decrease from top to bottom but not where the line crosses 50%')
    parser.add_argument('--combine32Returnliteral50', action='store_true', dest='combine32returnliteral50', help='use this flag to correct the motility ic50 (2-1-0 scoring) reporting if your ic50 fit reports the 50% decrease from top to bottom but not where the line crosses 50%')

    ##constraining parameters
    ###constraining hill slope
    parser.add_argument('--constrainMotHill', action='store_true', dest='constrainMotHill_bool', help='provide if and only if you want to constrain the hill slope for motility scoring (e.g., 3-2-1-0 scoring) and then use the --motHill argument to provide the constrained Hill Value')
    parser.add_argument('--motHill', action='store', dest='motHill', type=float, default=np.nan, required = False, help='the constant/constrained hill value that will be used rather than curve fitting attempting to optimize; required use when there are fewer than 3 experiments for a given drug/strain/stage treatment combo')
    parser.add_argument('--motnotDefHill', action='store', dest='motnotDefHill', type=float, default=np.nan, required=False, help='not default hill slope initial guess parameter for motility/IC50 analysis')
    parser.add_argument('--constrainMorHill', action='store_true', dest='constrainMorHill_bool', help='provide if and only if you want to constrain the hill slope for mortality scoring (e.g., 1-0 scoring) and then use the --morHill argument to provide the constrained Hill Value')
    parser.add_argument('--morHill', action='store', dest='morHill', type=float, default=np.nan, required = False, help='the constant/constrained hill value that will be used rather than curve fitting attempting to optimize; required use when there are fewer than 3 experiments for a given drug/strain/stage treatment combo')
    parser.add_argument('--mornotDefHill', action='store', dest='mornotDefHill', type=float, default=np.nan, required=False, help='not default Hill slope initial guess parameter for mortality/LC50 analysis')
    ### constraining top or bottom value
    parser.add_argument('--motnotDefTop', action='store', dest='motnotDefTop', type=float, default=np.nan, required=False, help='not default Top initial guess parameter for IC50 analysis with motility scoring. Default is the observed value in the dataset')
    parser.add_argument('--motnotDefBottom', action='store', dest='motnotDefBottom', type=float, default=np.nan, required=False, help='not default Bottom initial guess parameter for IC50 analysis with motility scoring. Default is the observed value in the dataset')
    parser.add_argument('--plotMotTopBottom', action='store_true', dest='plotMotTopBottom', help='use this flag if you want to plot with the true/observed top value rather than the fit top for motility scoring')
    parser.add_argument('--plotMotHundredZero', action='store_true', dest='plotMotHundredZero', help='use this flag if you want to plot with a top of 100 and a bottom of 0 no matter the true/observed top/bottom values or the curve fit values for motility scoring; do not use this flag with --plotMotTopBottom')

    parser.add_argument('--mornotDefTop', action='store', dest='mornotDefTop', type=float, default=np.nan, required=False, help='not default Top initial guess parameter for LC50 analysis with mortality scoring. Default is the observed value in the dataset')
    parser.add_argument('--mornotDefBottom', action='store', dest='mornotDefBottom', type=float, default=np.nan, required=False, help='not default Bottom initial guess parameter for LC50 analysis with mortality scoring. Default is the observed value in the dataset')
    parser.add_argument('--plotMorTopBottom', action='store_true', dest='plotMorTopBottom', help='use this flag if you want to plot with the true/observed top value rather than the fit top for mortality scoring')
    parser.add_argument('--plotMorHundredZero', action='store_true', dest='plotMorHundredZero', help='use this flag if you want to plot with a top of 100 and a bottom of 0 no matter the true/observed top/bottom values or the curve fit values for mortality scoring; do not use this flag with --plotMorTopBottom')

    ###guess ic50
    parser.add_argument('--motnotDefic50', action='store', dest='motnotDefic50', type=float, default=np.nan, required=False, help='not default ic50 initial guess parameter for motility/IC50 analysis')
    parser.add_argument('--mornotDefic50', action='store', dest='mornotDefic50', type=float, default=np.nan, required=False, help='not default ic50 initial guess parameter for mortality/LC50 analysis')


    ## spline orders if no fit is possible
    parser.add_argument('--motspline_k1', action='store', dest='motspline_k1', type=int, default=1, required=False, help='the order of the first part (small concentration) of the spline smoothing if there is no IC50 fit for the motility scoring (e.g., 3-2-1-0 scoring) analysis')
    parser.add_argument('--motspline_k2', action='store', dest='motspline_k2', type=int, default=1, required=False, help='the order of the second part (rest of concentrations) of the spline smoothing if there is no IC50 fit for the motility scoring (e.g., 3-2-1-0 scoring) analysis')
    parser.add_argument('--morspline_k1', action='store', dest='morspline_k1', type=int, default=1, required=False, help='the order of the first part (small concentration) of the spline smoothing if there is no LC50 fit for the mortality scoring (e.g., 3-2-1-0 scoring) analysis')
    parser.add_argument('--morspline_k2', action='store', dest='morspline_k2', type=int, default=1, required=False, help='the order of the second part (rest of concentrations) of the spline smoothing if there is no LC50 fit for the mortality scoring (e.g., 3-2-1-0 scoring) analysis')


    #arguments specific to probability model analyses if unequal worm numbers in well
    parser.add_argument('--random_seed', action='store', dest='random_seed', type=int, default=73)

    #arguments specific to stats analyses
    parser.add_argument('--stats_compare_to_control', action='store_true', dest='stats_compare_to_control', help='use this flag if and only if you want to run chi square comparisons to control; use --stats_concs and --stats_days to specify what days and concentrations are compared specifically; and make sure to include at least one or both of --stats_motility or --stats_mortality')
    parser.add_argument('--stats_make_every_possible_comparison', action='store_true', dest='stats_make_every_possible_comparison', help = 'use this flag if and only if you want to make every possible comparison in the chi-square analysis. Note this is NOT recommended as it increases the testing burden to reach significance due to multiple hypothesis testing correction')
    parser.add_argument('--stats_compare_to', action='store', nargs='+', dest='stats_compare_to', default=[0], help='all of the concentration indices you want to be the expecation/compared to. The default is just the lowest/control concentration. Add an integer for each additional comparison (e.g., --stats_compare_to 1 2 will also make comparisons with the 0.1 and 1 concentrations as the expectations); ignored if --stats_make_every_possible_comparison is provided')
    parser.add_argument('--stats_days', action='store', nargs='+', type=int, dest='stats_days', default=[1,4,7], help='which days to make the comparisons; ignored if --stats_make_every_possible_comparison is provided')
    parser.add_argument('--stats_concs', action='store', nargs='+', type=int, dest='stats_concs', default=[-1], help='the concentrations which will be compared to the expectation concentration. Default is the highest concentration; ignored if --stats_make_every_possible_comparison is provided. ')
    parser.add_argument('--stats_mortality', action='store_true', help='include if and only if you have included --stats_compare_to_control and want to test the differences in mortality; ignored if --stats_make_every_possible_comparison is provided')
    parser.add_argument('--stats_motility', action='store_true', help='include if and only if you have included --stats_compare_to_control and want to test the differences in motility, specifically inhibition profile (not highest score); ignored if --stats_make_every_possible_comparison is provided')
    return (parser)

def main():
    setup_formatting()
    parser_args = generate_parser().parse_args()
    analysis_instance = WormAnalysis(parser_args)
    analysis_instance.run_analysis()

class WormAnalysis():
    def __init__(self, parser_args):
        logfile='worm_analysis_{}_{}_{}.txt'.format(parser_args.drug, parser_args.stage, parser_args.strain)
        logging.basicConfig(filename=logfile, level=logging.INFO, filemode='w', format='%(name)s - %(levelname)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
        logging.info("command used to run analysis:\n {}".format(' '.join(sys.argv)))

        #metadata
        self.drug = parser_args.drug
        self.strain = parser_args.strain
        self.stage = parser_args.stage
        self.concUnits = worm_analysis_utils.find_concUnits(parser_args.concUnits)
        self.molarmass = parser_args.molarmass
        self.density = parser_args.density

        #data of interest
        self.C_day = parser_args.C_day
        self.x0_value = parser_args.x0_val
        self.motreturnliteral50 = parser_args.motreturnliteral50
        self.morreturnliteral50 = parser_args.morreturnliteral50
        self.combine32returnliteral50 = parser_args.combine32returnliteral50
        self.representative = parser_args.representative
        self.num_replicates = parser_args.num_replicates
        self.nscore_insystem = 4

        #analyses to do
        self.reportNum = parser_args.reportNum
        self.linePlots = parser_args.linePlots
        self.linePlotsMot = parser_args.linePlotsMot
        self.linePlotsMor = parser_args.linePlotsMor
        self.include_single_exp_plots = parser_args.isep
        self.include_sem = parser_args.isem
        self.include_stdev = parser_args.istdev
        self.time50 = parser_args.time50
        self.time50Mot = parser_args.plotIT50
        self.time50Mor = parser_args.plotLT50
        self.conc50 = parser_args.conc50
        self.conc50Mot = parser_args.plotIC50
        self.conc50Mor = parser_args.plotLC50
        self.runNo3 = parser_args.runNo3

        #curve fitting and plotting
        self.constrainNo3Hill, self.constrainMorHill, self.constrainMotHill = parser_args.constrainNo3Hill_bool, parser_args.constrainMorHill_bool,  parser_args.constrainMotHill_bool
        self.no3Hill, self.morHill, self.motHill = parser_args.no3Hill, parser_args.morHill, parser_args.motHill

        self.motspline_k1, self.motspline_k2 = parser_args.motspline_k1, parser_args.motspline_k2
        self.morspline_k1, self.morspline_k2 = parser_args.morspline_k1, parser_args.morspline_k2
        self.no3spline_k1, self.no3spline_k2 = parser_args.no3spline_k1, parser_args.no3spline_k2

        self.motnotDefTop, self.motnotDefBottom, self.motnotDefic50, self.motnotDefHill  = parser_args.motnotDefTop, parser_args.motnotDefBottom, parser_args.motnotDefic50, parser_args.motnotDefHill
        self.mornotDefTop, self.mornotDefBottom, self.mornotDefic50, self.mornotDefHill = parser_args.mornotDefTop, parser_args.mornotDefBottom, parser_args.mornotDefic50, parser_args.mornotDefHill
        self.no3notDefHill, self.no3notDefTop, self.no3notDefBottom, self.no3notDefic50 = parser_args.no3notDefHill, parser_args.no3notDefTop, parser_args.no3notDefBottom, parser_args.no3notDefic50


        self.plotMotTopBottom_bool, self.plotMorTopBottom_bool, self.plotno3TopBottom_bool = parser_args.plotMotTopBottom, parser_args.plotMorTopBottom, parser_args.plotno3TopBottom
        self.plotMotHundredZero, self.plotMorHundredZero, self.plotno3HundredZero = parser_args.plotMotHundredZero, parser_args.plotMorHundredZero, parser_args.plotno3HundredZero



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

    def make_tables(self):
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


    def load_data(self, toAssess, expNames):
        self.num_experiments = len(toAssess)
        self.consistent_num = np.tile(np.nan, self.num_experiments)
        self.expNames = expNames

        scores_arrs_conc = {}
        scores_arrs_well = {}
        for i, file in enumerate(toAssess):
            df = pd.read_csv(file, sep='\t', header=0)
            if i == 0:
                self.uniq_conc = np.unique(df.loc[:,'Concentration'])
                self.mM = worm_analysis_utils.find_mM(self.uniq_conc, self.concUnits, self.molarmass, self.density)
                self.num_concentrations = np.shape(self.uniq_conc)[0]
                self.num_days = np.shape(np.unique(df.loc[:,'Day']))[0]
                assert 0 not in np.unique(np.array(df.loc[:,'Day']).astype(np.int32)), "Remove Day 0 info from input file"

                self.make_tables()

                '''set some dictionaries for mapping later'''
                self.conc_to_index = {float(c):k for k,c in enumerate(self.uniq_conc)}
                self.index_to_conc = {k:float(c) for k,c in enumerate(self.uniq_conc)}

                self.well_index_to_conc = {w-1: float(c) for w,c in zip(np.array(df.loc[:,'Well']).astype(np.int32)[0:self.num_concentrations*self.num_replicates], np.array(df.loc[:,'Concentration'])[0:self.num_concentrations*self.num_replicates])}
                self.conc_to_well_index = {float(c): w-1 for w,c in zip(np.array(df.loc[:,'Well']).astype(np.int32)[0:self.num_concentrations*self.num_replicates], np.array(df.loc[:,'Concentration'])[0:self.num_concentrations*self.num_replicates])}


            scores_arrs_conc[i] = self.find_scores_fromdf(df, conc=True)
            scores_arrs_well[i] = self.find_scores_fromdf(df, well=True)
            consistent_num = np.nansum(df.loc[:,"numTotal_equal"] == True) == self.num_concentrations*self.num_replicates
            self.consistent_num[i] = consistent_num

        self.scores3_by_conc = np.stack(([scores_arrs_conc[i] for i in np.arange(len(toAssess))]), axis=-1)
        self.scores3_by_well = np.stack(([scores_arrs_well[i] for i in np.arange(len(toAssess))]), axis=-1)


        logging.info("Loaded all the data")
        if not np.all(self.consistent_num):
            logging.warning('wells in {} do not have equal numbers of worms in a given well across the length of the experiment'.format(np.array(expNames)[np.where(self.consistent_num == False)[0]]))
            worm_analysis_utils.plot_wormnums(self.scores3_by_well, self.expNames, self.num_days)

    def drive_linePlots(self, figname_base=''):
        if self.linePlots:
            if self.linePlotsMot:
                if self.include_single_exp_plots:
                    for (i, exp) in enumerate(self.expNames):
                        worm_analysis_scorelines.plotLineCurves(self.motility_index_scores_by_conc[:,:,i], 'isep_motility_{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, exp, figname_base), self.expNames[i], "Motility Index Score", 0, self.nscore_insystem-1, 1, np.arange(self.num_days + 1), self.uniq_conc, self.concUnits, self.mM)
                to_plot_mot = worm_analysis_scorelines.get_avg_across_exp(self.motility_index_scores_by_well, self.num_concentrations, self.num_replicates, self.num_experiments, self.num_days)
                worm_analysis_scorelines.plotLineCurves(to_plot_mot,'average_motility_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, figname_base), r"\textbf{%s on %s %s}" %(self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$", "Motility Index Score", 0, self.nscore_insystem-1, 1, np.arange(self.num_days + 1), self.uniq_conc, self.concUnits, self.mM)
            if self.linePlotsMor:
                if self.include_single_exp_plots:
                    for i, exp in enumerate(self.expNames):
                        worm_analysis_scorelines.plotLineCurves(self.mortality_scores_by_conc[:,:,i], 'isep_mortality_{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, exp, figname_base), exp, "\% Alive", 0, 100, 25, np.arange(self.num_days + 1), self.uniq_conc, self.concUnits, self.mM)
                to_plot_mor = worm_analysis_scorelines.get_avg_across_exp(self.mortality_scores_by_well, self.num_concentrations, self.num_replicates, self.num_experiments, self.num_days)
                worm_analysis_scorelines.plotLineCurves(to_plot_mor,'average_mortality_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, figname_base), r"\textbf{%s on %s %s}" %(self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$" , "\% Alive", 0, 100, 25, np.arange(self.num_days + 1), self.uniq_conc, self.concUnits, self.mM)

    def drive_survivalTimePlots(self, figname_base='', no3=False):
        if self.time50:
            if self.representative  == 0:
                logging.error("Must provide a representative experiment for Survival Analysis using --representative argument")
                exit(1)
            logging.info('Beginning Survival Analysis for {} as the representative experiment'.format(self.expNames[self.representative-1]))
            if self.time50Mot:
                inhibited = np.full((self.num_concentrations, self.num_days+1), 100, dtype=np.float64)
                inhibited, logging_value, IT50 = worm_analysis_it.findSurvivability(inhibited, self.scores3_by_conc, self.representative-1, self.expNames[self.representative-1], self.nscore_insystem, self.consistent_num[self.representative-1], self.num_concentrations, self.num_days, motility=True, no3=no3)
                self.df_tablec[logging_value] = IT50
                logging.info("Added the IT50 values to the table which will be printed later. Column name is 'IT50'")
                if no3:
                    inhibited_og = np.full((self.num_concentrations, self.num_days+1), 100, dtype=np.float64)
                    inhibited_og, logging_value_og, IT50_og = worm_analysis_it.findSurvivability(inhibited_og, self.og_scores3_by_conc, self.representative-1, self.expNames[self.representative-1], self.nscore_insystem+1, self.consistent_num[self.representative-1], self.num_concentrations, self.num_days, motility=True)
            if self.time50Mor:
                mortality = np.full((self.num_concentrations, self.num_days+1), 100, dtype=np.float64)
                mortality, logging_value, LT50 = worm_analysis_it.findSurvivability(mortality, self.scores3_by_conc, self.representative-1, self.expNames[self.representative-1], self.nscore_insystem, self.consistent_num[self.representative-1], self.num_concentrations, self.num_days, mortality=True)
                self.df_tablec[logging_value] = LT50
                logging.info("Added the LT50 values to the table which will be printed later. Column name is 'LT50'")
            if self.time50Mot and self.time50Mor:
                figname_base_spec = '{}_IT50_LT50_' + '{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, self.expNames[self.representative-1], figname_base)
                if no3:
                    worm_analysis_it.plotSurvivalTime(inhibited, mortality, figname_base_spec, self.uniq_conc, self.concUnits, self.drug, self.stage, self.strain, np.arange(self.num_days + 1), no3=no3, inhibited_og = inhibited_og)
                else:
                    worm_analysis_it.plotSurvivalTime(inhibited, mortality, figname_base_spec, self.uniq_conc, self.concUnits, self.drug, self.stage, self.strain, np.arange(self.num_days + 1))
            else:
                if self.time50Mot:
                    figname_base_spec = '{}_IT50_' + '{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, self.expNames[self.representative-1], figname_base)
                    worm_analysis_it.plotSurvivalTime(inhibited, np.nan, figname_base_spec, self.uniq_conc, self.concUnits, self.drug, self.stage, self.strain, np.arange(self.num_days + 1), plot_mortality=False)
                if self.time50Mor:
                    logging.warning('Why do you want to plot only the LT50? I suggest plotting IT50 and LT50 together. But here you go')
                    figname_base_spec = '{}_LT50_' + '{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, self.expNames[self.representative-1], figname_base)
                    worm_analysis_it.plotSurvivalTime(np.nan, mortality, figname_base_spec, self.uniq_conc, self.concUnits, self.drug, self.stage, self.strain, np.arange(self.num_days +1), plot_motility=False)
        logging.info('Completed Survival Analysis')


    def transform_X(self, X):
        '''set [0] to x0_val'''
        bool_0 = X == 0
        X[bool_0] += self.x0_value
        '''return log10 of concentrations'''
        return np.log10(X)

    def inhibitorResponse_equation(self, X, top, bottom, ic50, hillslope):
        '''
        Y=Bottom + (Top-Bottom)/(1+10^((LogIC50-X)*HillSlope))
        where X is log10([]) with x=0 -> x=1e-6 (i.e. call self.transform_X(X))
        '''
        exponent = (np.log10(ic50)- self.transform_X(X))*hillslope
        equation = bottom + (top-bottom)/(1+(10**exponent))
        return equation

    def inhibitorResponse_equation_constrainHillMor(self, X, top, bottom, ic50):
        '''
        Y=Bottom + (Top-Bottom)/(1+10^((LogIC50-X)*HillSlope))
        where X is log10([]) with x=0 -> x=1e-6 (i.e. call self.transform_X(X))
        '''
        exponent = (np.log10(ic50)- self.transform_X(X))*self.morHill
        equation = bottom + (top-bottom)/(1+(10**exponent))
        return equation

    def inhibitorResponse_equation_constrainHillMot(self, X, top, bottom, ic50):
        '''
        Y=Bottom + (Top-Bottom)/(1+10^((LogIC50-X)*HillSlope))
        where X is log10([]) with x=0 -> x=1e-6 (i.e. call self.transform_X(X))
        '''
        exponent = (np.log10(ic50)- self.transform_X(X))*self.motHill
        equation = bottom + (top-bottom)/(1+(10**exponent))
        return equation

    def inhibitorResponse_equation_constrainHillNo3(self, X, top, bottom, ic50):
        '''
        Y=Bottom + (Top-Bottom)/(1+10^((LogIC50-X)*HillSlope))
        where X is log10([]) with x=0 -> x=1e-6 (i.e. call self.transform_X(X))
        '''
        exponent = (np.log10(ic50)- self.transform_X(X))*self.no3Hill
        equation = bottom + (top-bottom)/(1+(10**exponent))
        return equation

    def drive_inhibitoryConcentrationPlots(self, figname_base='', no3=False):
        if self.conc50:
            conc_X = np.tile(np.array([self.well_index_to_conc[x] for x in np.arange(self.num_concentrations*self.num_replicates)]).reshape(-1, 1), (1, self.num_experiments))
            if self.conc50Mot:
                #find percent inhibited each concentration for that day of analysis (self.C_day), averages, and standard error around the means
                uninhibitedMot, avgMot, semMot = worm_analysis_ic.find_avg_sem(self.scores3_by_well, self.C_day, self.nscore_insystem, self.num_concentrations, self.num_replicates, self.num_experiments, self.conc_to_index, self.well_index_to_conc, motility = True)
                #set initial guesses for top, bottom, hill slope, and ic50
                if no3:
                    mot_initial_guesses = worm_analysis_ic.set_guess_params(self.num_experiments, avgMot, self.no3notDefHill, self.no3notDefTop, self.no3notDefBottom, self.no3notDefic50, constrainHill = self.constrainNo3Hill, motility=True)
                    #curve fit based on initial guess params
                    logging.info('Running Levenberg-Marquardt Algorithm Scipy Curve Fitting for no3 Motility (2-1-0) scoring using a max number of function evaluations of {}. Initial values are the following.\nINITIAL Top:\t{}\nINITIAL Bottom:\t{}\nINITIAL IC50:\t{}\nINITIAL HillSlope:\t{}'.format(int(1e6), mot_initial_guesses["P0_top"], mot_initial_guesses["P0_bottom"], mot_initial_guesses["P0_ic50"], mot_initial_guesses["P0_hill"]))
                    if self.num_experiments < 3 or self.constrainNo3Hill:
                        popt, popc = curve_fit(self.inhibitorResponse_equation_constrainHillNo3, conc_X.flatten(), uninhibitedMot.flatten(), p0=mot_initial_guesses["P0"],method='lm', maxfev=int(1e6))
                    else:
                        popt, popc = curve_fit(self.inhibitorResponse_equation, conc_X.flatten(), uninhibitedMot.flatten(), p0=mot_initial_guesses["P0"], method='lm', maxfev=int(1e6))
                    logging.info("Returned fit for no3 Motility (2-1-0) scoring (top, bottom, ic50, hill (if number of experiments >3)): {}".format(popt))
                    nofit_bool = worm_analysis_ic.evaluate_no_fit(popt[0], popt[1], popt[2], np.sort(self.uniq_conc)[-1])
                    #plot curve fit
                    mot_splineic50 = worm_analysis_ic.plotIC(r'$\mathrm{IC_{50}}$' + ' {} on {} {} Day {}'.format(self.drug, self.stage, self.strain, self.C_day), 'IC50_{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, self.C_day, figname_base), avgMot, semMot, self.uniq_conc, self.x0_value, self.concUnits, self.no3spline_k1, self.no3spline_k2, popt, constrainedHill = self.no3Hill, use100_0=self.plotno3HundredZero, useobserved=self.plotno3TopBottom_bool, fit_possible = not nofit_bool, returnliteral50 = self.combine32returnliteral50)
                    #record value in table
                    if nofit_bool or self.combine32returnliteral50:
                        self.df_tabled = worm_analysis_ic.set_to_log_value(self.df_tabled, self.C_day, popt[1], mot_splineic50, np.sort(self.uniq_conc)[-1], np.sort(self.uniq_conc)[1], motility=True, no3=no3)
                    else:
                        self.df_tabled = worm_analysis_ic.set_to_log_value(self.df_tabled, self.C_day, popt[1], popt[2], np.sort(self.uniq_conc)[-1], np.sort(self.uniq_conc)[1], motility=True, no3=no3)

                else:
                    mot_initial_guesses = worm_analysis_ic.set_guess_params(self.num_experiments, avgMot, self.motnotDefHill, self.motnotDefTop, self.motnotDefBottom, self.motnotDefic50, constrainHill = self.constrainMotHill, motility=True)
                    #curve fit based on initial guess params
                    logging.info('Running Levenberg-Marquardt Algorithm Scipy Curve Fitting for Motility (3-2-1-0) scoring using a max number of function evaluations of {}. Initial values are the following.\nINITIAL Top:\t{}\nINITIAL Bottom:\t{}\nINITIAL IC50:\t{}\nINITIAL HillSlope:\t{}'.format(int(1e6), mot_initial_guesses["P0_top"], mot_initial_guesses["P0_bottom"], mot_initial_guesses["P0_ic50"], mot_initial_guesses["P0_hill"]))
                    if self.num_experiments < 3 or self.constrainMotHill:
                        popt, popc = curve_fit(self.inhibitorResponse_equation_constrainHillMot, conc_X.flatten(), uninhibitedMot.flatten(), p0=mot_initial_guesses["P0"], method='lm', maxfev=int(1e6))
                    else:
                        popt, popc = curve_fit(self.inhibitorResponse_equation, conc_X.flatten(), uninhibitedMot.flatten(), p0=mot_initial_guesses["P0"], method='lm', maxfev=int(1e6))
                    logging.info("Returned fit for Motility (3-2-1-0) scoring (top, bottom, ic50, hill (if number of experiments >3)): {}".format(popt))
                    nofit_bool = worm_analysis_ic.evaluate_no_fit(popt[0], popt[1], popt[2], np.sort(self.uniq_conc)[-1])
                    #plot curve fit
                    mot_splineic50 = worm_analysis_ic.plotIC(r'$\mathrm{IC_{50}}$' + ' {} on {} {} Day {}'.format(self.drug, self.stage, self.strain, self.C_day), 'IC50_{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, self.C_day, figname_base), avgMot, semMot, self.uniq_conc, self.x0_value, self.concUnits, self.motspline_k1, self.motspline_k2, popt, constrainedHill = self.motHill, use100_0=self.plotMotHundredZero, useobserved=self.plotMotTopBottom_bool, fit_possible = not nofit_bool, returnliteral50 = self.motreturnliteral50)
                    #record value in table
                    if nofit_bool or self.motreturnliteral50:
                        self.df_tabled = worm_analysis_ic.set_to_log_value(self.df_tabled, self.C_day, popt[1], mot_splineic50, np.sort(self.uniq_conc)[-1], np.sort(self.uniq_conc)[1], motility=True)
                    else:
                        self.df_tabled = worm_analysis_ic.set_to_log_value(self.df_tabled, self.C_day, popt[1], popt[2], np.sort(self.uniq_conc)[-1], np.sort(self.uniq_conc)[1], motility=True)
            if self.conc50Mor:
                uninhibitedMor, avgMor, semMor = worm_analysis_ic.find_avg_sem(self.scores3_by_well, self.C_day, self.nscore_insystem, self.num_concentrations, self.num_replicates, self.num_experiments, self.conc_to_index, self.well_index_to_conc, mortality = True)
                mor_initial_guesses = worm_analysis_ic.set_guess_params(self.num_experiments, avgMor, self.mornotDefHill, self.mornotDefTop, self.mornotDefBottom, self.mornotDefic50, constrainHill = self.constrainMorHill, mortality=True)
                logging.info('Running Levenberg-Marquardt Algorithm Scipy Curve Fitting for Mortality (1-0) scoring using a max number of function evaluations of {}. Initial values are the following.\nINITIAL Top:\t{}\nINITIAL Bottom:\t{}\nINITIAL IC50:\t{}\nINITIAL HillSlope:\t{}'.format(int(1e6), mor_initial_guesses["P0_top"], mor_initial_guesses["P0_bottom"], mor_initial_guesses["P0_ic50"], mor_initial_guesses["P0_hill"]))
                if self.num_experiments < 3 or self.constrainMorHill:
                    popt, popc = curve_fit(self.inhibitorResponse_equation_constrainHillMor, conc_X.flatten(), uninhibitedMor.flatten(), p0=mor_initial_guesses["P0"], method='lm', maxfev=int(1e6))
                else:
                    popt, popc = curve_fit(self.inhibitorResponse_equation, conc_X.flatten(), uninhibitedMor.flatten(), p0=mor_initial_guesses["P0"], method='lm', maxfev=int(1e6))
                logging.info("Returned fit for Mortality (1-0) scoring (top, bottom, lc50, hill (if number of experiments >3)): {}".format(popt))
                nofit_bool = worm_analysis_ic.evaluate_no_fit(popt[0], popt[1], popt[2], np.sort(self.uniq_conc)[-1])
                #plot curve fit
                mor_splineic50 = worm_analysis_ic.plotIC(r'$\mathrm{LC_{50}}$' + ' {} on {} {} Day {}'.format(self.drug, self.stage, self.strain, self.C_day), 'LC50_{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, self.C_day, figname_base), avgMor, semMor, self.uniq_conc, self.x0_value, self.concUnits, self.morspline_k1, self.morspline_k2, popt, constrainedHill = self.morHill, use100_0=self.plotMorHundredZero, useobserved=self.plotMorTopBottom_bool, fit_possible = not nofit_bool, returnliteral50 = self.morreturnliteral50)
                #record value in table
                if nofit_bool or self.morreturnliteral50:
                    self.df_tabled = worm_analysis_ic.set_to_log_value(self.df_tabled, self.C_day, popt[1], mor_splineic50, np.sort(self.uniq_conc)[-1], np.sort(self.uniq_conc)[1], mortality=True)
                else:
                    self.df_tabled = worm_analysis_ic.set_to_log_value(self.df_tabled, self.C_day, popt[1], popt[2], np.sort(self.uniq_conc)[-1], np.sort(self.uniq_conc)[1], mortality=True)
            logging.info('Completed Non-linear Regression for Inhibition Response Analysis')


    def drive_tablereports(self, figname_base = ''):
        worm_analysis_report.reportTable(self.df_tablec, self.df_tabled, self.representative-1, self.expNames, self.C_day, self.drug, self.stage, self.strain, self.reportNum, self.time50Mot, self.time50Mor, self.conc50Mot, self.conc50Mor, figname_base=figname_base)


    def reset_to_no3_system(self):
        self.nscore_insystem = 3

        self.og_scores3_by_conc, self.og_scores3_by_well = self.scores3_by_conc.copy(), self.scores3_by_well.copy()
        self.og_motility_index_scores_by_conc, self.og_motility_index_scores_by_well = self.motility_index_scores_by_conc.copy(), self.motility_index_scores_by_well.copy()
        self.scores3_by_conc, self.scores3_by_well = worm_analysis_no3.transform_raw_to_no3(self.og_scores3_by_conc, self.nscore_insystem, self.num_days, self.num_experiments), worm_analysis_no3.transform_raw_to_no3(self.og_scores3_by_well, self.nscore_insystem, self.num_days, self.num_experiments)
        self.motility_index_scores_by_conc, self.motility_index_scores_by_well = worm_analysis_utils.find_scores_driver(self.scores3_by_conc, self.scores3_by_well, self.num_concentrations, self.num_days, self.num_experiments, self.num_replicates, self.nscore_insystem, motility = True)
        logging.info('Transformed data by combining all scores of 2&3 for each day, concentration, and experiment & found the new motility index scores')


        self.oglinePlotsMor, self.ogtime50Mor, self.ogconc50Mor = self.linePlotsMor, self.time50Mor, self.conc50Mor
        self.linePlotsMor, self.conc50Mor = False, False

        self.df_tabled[r'$\mathrm{\textbf{IC50}\;(combined\;3\;\&\;2)}$'] = np.tile('NC', self.df_tabled.shape[0])
        #self.make_tables()

    def driveStats(self):
        if self.stats_compare_to_control:
            pval_dfs = {}
            if self.stats_make_every_possible_comparison:
                for conc in self.uniq_conc:
                    compare_to_index = self.conc_to_index[conc]
                    compare_with_concs = self.uniq_conc[self.uniq_conc != conc]
                    df_table_i, df_table_m = pd.DataFrame(index=compare_with_concs), pd.DataFrame(index=compare_with_concs)
                    for df in [df_table_i, df_table_m]:
                        df.index.name = "Concentration ({})".format(self.concUnits)
                    for compare_to_day in np.arange(1, self.num_days + 1):
                        pvals_i = np.array([])
                        pvals_m = np.array([])
                        expected_inh = worm_analysis_stats.make_inh_arr(self.scores3_by_conc, compare_to_index, compare_to_day-1, self.nscore_insystem)
                        expected_mor = worm_analysis_stats.make_mor_arr(self.scores3_by_conc, compare_to_index, compare_to_day-1)
                        for compare_with_conc in compare_with_concs:
                            compare_with_conc_index = self.conc_to_index[compare_with_conc]
                            compare_with_inh = worm_analysis_stats.make_inh_arr(self.scores3_by_conc, compare_with_conc_index, compare_to_day-1, self.nscore_insystem)
                            compare_with_mor = worm_analysis_stats.make_mor_arr(self.scores3_by_conc, compare_with_conc_index, compare_to_day-1)
                            raw_pval_inh = worm_analysis_stats.chisquare_it(compare_with_inh, expected_inh)
                            raw_pval_mor = worm_analysis_stats.chisquare_it(compare_with_mor, expected_mor)
                            pvals_i = np.hstack((pvals_i, raw_pval_inh))
                            pvals_m = np.hstack((pvals_m, raw_pval_mor))
                        df_table_i["Day {}".format(compare_to_day)] = pvals_i
                        df_table_m["Day {}".format(compare_to_day)] = pvals_m
                    pval_dfs[conc] = {"inh": df_table_i, "mor": df_table_m}
                pval_dfs, num_tests = worm_analysis_stats.multiple_hypothesis_testing_correction(df_dict)
            else:
                #what about specified testing only?
                filler = 0
        plot_it() #need to make this function

    def reset_to_og_score_system(self):
        self.nscore_insystem = 4

        self.scores3_by_conc, self.scores3_by_well = self.og_scores3_by_conc.copy(), self.og_scores3_by_well.copy()
        self.motility_index_scores_by_conc, self.motility_index_scores_by_well = self.og_motility_index_scores_by_conc.copy(), self.og_motility_index_scores_by_well.copy()

        self.linePlotsMor, self.time50Mor, self.conc50Mor = self.oglinePlotsMor, self.ogtime50Mor, self.ogconc50Mor

        self.make_tables()

    def drive_no3analysis(self):
        if self.runNo3:
            self.reset_to_no3_system()
            logging.info('restarting the analyses with the combine 3-2 scoring system and new motility index score')
            self.drive_linePlots(figname_base = "_combine3-2_analysis")
            self.drive_survivalTimePlots(figname_base = "_combine3-2_analysis", no3=self.runNo3)
            self.drive_inhibitoryConcentrationPlots(figname_base = "_combine3-2_analysis", no3=self.runNo3)
            ## drive stats
            self.drive_tablereports(figname_base = '_combine3-2_analysis')
            self.reset_to_og_score_system()

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
        self.drive_inhibitoryConcentrationPlots()

        # Drive stats and table reports
        ## drive stats
        self.drive_tablereports()

        # Drive no3 score analyses
        self.drive_no3analysis()

        # Drive unequal worms correction

        ## repeat line plots, surivival analysis, inhibitory concentration analysis, no3 analyses also

main()
