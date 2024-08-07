#!/usr/bin/env python3

'''
Usage example:
./screening_analysis.py --toAssess test_files/kw9_export.txt test_files/kw11_export.txt test_files/kw12_export.txt --expNames kw9 kw11 kw12 --strain N2 --lifestage L4 --drug ALB --concUnits 0 --representative 3 --molarmass 265.333
'''

import numpy as np
import argparse as ap
from scipy.optimize import curve_fit
from scipy.stats import sem, chisquare
import logging
import pandas as pd
import datetime
import matplotlib.pyplot as plt
import sys
from matplotlib import rc, rcParams
from scipy.interpolate import make_interp_spline, BSpline
from screening_analysis_no3 import WormAnalysis_no3
from screening_analysis_prob_model import WormAnalysis_probModel

rc('axes', linewidth=2)
params = {'font.sans-serif': 'Helvetica',
          'font.size': 12,
          'font.weight': 'bold',
          'legend.frameon': False,
          'legend.labelspacing': 1,
          'text.usetex': True,
          'text.latex.preamble': [r'\usepackage{siunitx}',
                                  r'\sisetup{detect-all}',
                                  r'\usepackage{sansmath}',
                                  r'\sansmath']}
rcParams.update(params)

def main():
    parser = generate_parser()
    args = parser.parse_args()
    analysis_instance_main = WormAnalysis(args.toAssess, args.expNames, args.drug, args.strain, args.stage, args.concUnits, args.molarmass, args.density, args.reportNum)
    analysis_structure(args, analysis_instance_main)

    if args.runNo3:
        logging.info('now beginning the analysis where scores 2&3 are combined.')
        analysis_instance_no3 = WormAnalysis_no3(analysis_instance_main.consistent_num, analysis_instance_main.drug, analysis_instance_main.strain, analysis_instance_main.stage, analysis_instance_main.uniq_conc, analysis_instance_main.concUnits, analysis_instance_main.concUnits_dict, analysis_instance_main.conc_index, analysis_instance_main.well_index_to_conc, analysis_instance_main.conc_colors_lo_to_hi, analysis_instance_main.conc_markers_lo_to_hi, analysis_instance_main.conc_marker_outline_lo_to_hi, analysis_instance_main.C_day, analysis_instance_main.x0_val, analysis_instance_main.mM, analysis_instance_main.num_days, analysis_instance_main.num_experiments, analysis_instance_main.num_concentrations, analysis_instance_main.scores3_by_well, analysis_instance_main.scores3_by_conc)
        analysis_instance_no3.run(args.plotLine3, args.plotIT50, args.plotIC50, args.isep, args.expNames, args.rep, args.C_day, args.x0_val, args.hill2, args.constrain2Hill_bool, args.spline_k1, args.spline_k2)


    if not analysis_instance_main.consistent_num:
        print("This prob model part is almost ready, but not verified quite yet")
        #logging.info('now beginning the analysis using a probability model since an inconsistent number of worms were observed in atleast one well across the duration of the experiment.')
        #analysis_instance_probModel = WormAnalysis_probModel(analysis_instance_main.drug, analysis_instance_main.strain, analysis_instance_main.stage, analysis_instance_main.num_days, analysis_instance_main.num_experiments, analysis_instance_main.num_concentrations, analysis_instance_main.scores3_by_well, analysis_instance_main.scores3_by_conc, analysis_instance_main.get_num_treated(), analysis_instance_main.total_num_in_well, analysis_instance_main.conc_index, analysis_instance_main.well_index_to_conc, analysis_instance_main.random_seed)
        # '''return new/full generated data'''
        #generated_scores3_by_well, generated_scores3_by_conc = analysis_instance_probModel.run()
        # '''something to change the self values here in analysis_instance_main'''
        #analysis_instance_main.reset_self_arrays(generated_scores3_by_well, generated_scores3_by_conc)
        #analysis_structure(args, analysis_instance_main, figname_base='_withProbModel')

        if args.runNo3: #want to fix the figname base -- do this
            logging.info('now beginning the analysis where scores 2&3 are combined.')
            analysis_instance_no3 = WormAnalysis_no3(analysis_instance_main.consistent_num, analysis_instance_main.drug, analysis_instance_main.strain, analysis_instance_main.stage, analysis_instance_main.uniq_conc, analysis_instance_main.concUnits, analysis_instance_main.concUnits_dict, analysis_instance_main.conc_index, analysis_instance_main.well_index_to_conc, analysis_instance_main.conc_colors_lo_to_hi, analysis_instance_main.conc_markers_lo_to_hi, analysis_instance_main.conc_marker_outline_lo_to_hi, analysis_instance_main.C_day, analysis_instance_main.x0_val, analysis_instance_main.mM, analysis_instance_main.num_days, analysis_instance_main.num_experiments, analysis_instance_main.num_concentrations, analysis_instance_main.scores3_by_well, analysis_instance_main.scores3_by_conc)
            analysis_instance_no3.run(args.plotLine3, args.plotIT50, args.plotIC50, args.isep, args.expNames, args.rep, args.C_day, args.x0_val, args.hill2, args.constrain2Hill_bool, args.spline_k1, args.spline_k2)



def analysis_structure(args, analysis_instance, figname_base=''):
    if args.plotLine3 or args.plotLine1:
        analysis_instance.driveLinePlots(args.plotLine3, args.plotLine1, args.isep, args.expNames, figname_base=figname_base)
    if args.plotIT50 or args.plotLT50:
        analysis_instance.driveSurvivalTimePlots(args.plotIT50, args.plotLT50, args.rep, args.expNames, figname_base=figname_base)
    if args.plotIC50 or args.plotLC50:
        if args.notDefaultC50:
            analysis_instance.driveIC(args.plotIC50, args.plotLC50, args.C_day, args.x0_val, args.notDefHill1, args.notDefHill3, args.spline_k1, args.spline_k2, args.fitTopPlot_bool, args.constrain3Hill_bool, args.constrain1Hill_bool, default=False, not_default_3 = {'top':args.notDefTop3,'bottom':args.notDefBot3,'ic50':args.notDefic50,'hill':args.notDefHill3}, not_default_1 = {'top':args.notDefTop1,'bottom':args.notDefBot1,'ic50':args.notDeflc50,'hill':args.notDefHill1}, figname_base=figname_base)
        else:
            analysis_instance.driveIC(args.plotIC50, args.plotLC50, args.C_day, args.x0_val, args.hill1, args.hill3, args.spline_k1, args.spline_k2, args.fitTopPlot_bool, args.constrain3Hill_bool, args.constrain1Hill_bool, figname_base=figname_base)
    if args.stats_compare_to_control:
        try:
            assert args.stats_inhibited or args.stats_mortality
        except AssertionError as err:
            logging.exception("Assertion failed: trying to run statistics but don't provide what basis you want to compare (--stats_inhibited and/or --stats_mortality); see documentationn for more information")
            raise err
        analysis_instance.driveStats(args.stats_compare_to, args.stats_days, args.stats_inhibited, args.stats_mortality)
    analysis_instance.reportTable(args.expNames[args.rep-1], args.reportNum, args.plotIT50, args.plotLT50, args.plotIC50, args.plotLC50, figname_base=figname_base)


def generate_parser():
    parser = ap.ArgumentParser(description='C elegans analysis')
    parser.add_argument('--toAssess', action='store', dest='toAssess', nargs='+', type=str, required=True, help='the exported files to assess')
    parser.add_argument('--expNames', action='store', dest='expNames', nargs='+', type=str, required=True, help='the names of the experiments passed in toAssess (same order); will be used to annotate plots and tables')
    parser.add_argument('--strain', action='store', dest='strain', type=str, required=True, help='the strain which was treated in the assays (ex. N2 or Hawaii); will be used to annotate plots and tables')
    parser.add_argument('--lifestage', action='store', dest='stage', type=str, required=True, help='the lifestage which was treated in the assays (ex. L1 or L4); will be used to annotate plots and tables')
    parser.add_argument('--drug', action='store', dest='drug', type=str, required=True, help='the drug used in treatment in the assays (ex. ALB, PYR, IVM, or NTZ, etc); will be used to annotate plots and tables; for now this only accepts a single value')
    parser.add_argument('--molarmass', action='store', dest='molarmass', type=float, required=True, help='the molar mass of the drug used in treatment in the assays. Must be given in the g/mol units')
    parser.add_argument('--densityLiq', action='store', dest='density', type=float, default=1.0, help='the density of the liquid; must be entered in g/mL; default is that of water (1 g/mL)')
    parser.add_argument('--concUnits', action='store', dest='concUnits', type=int, required=True, help='use to specify the units of the concentration. 0 is ug/mL, 1 is uM, 2 is M, 3 is mM, 4 is nM, 5 is ng/mL, 6 is mg/mL, 7 is mg/kg, 8 is mg/g')
    parser.add_argument('--no_plotLine3', action='store_false', dest='plotLine3', help='add this argument to skip plotting the daily motility response by concentration')
    parser.add_argument('--no_plotLine1', action='store_false', dest='plotLine1', help='add this argument to skip plotting the daily lethality response by concentration')
    parser.add_argument('--include_single_exp_plots', action='store_true', dest='isep', help='when plotting the daily motility or lethalty response by concentration, add this argument to plot single experiments as well as the average')
    parser.add_argument('--no_plotIC50', action='store_false', dest='plotIC50', help='add this argument to skip plotting the IC50 (3-2-1-0 scoring)')
    parser.add_argument('--no_plotLC50', action='store_false', dest='plotLC50', help='add this argument to skip plotting the LC50 (1-0 scoring)')
    parser.add_argument('--C_day', action='store', dest='C_day', type=int, default=4, help='the day (index from 1) to assess for inhibitory and/or lethal concentration')
    parser.add_argument('--x0_value', action='store', dest='x0_val', type=float, default=1e-6, help='value to replace the x=0 [] with when transforming x')
    parser.add_argument('--no_plotIT50', action='store_false', dest='plotIT50', help='add this argument to skip plotting the IT50 (3-2-1-0 scoring)')
    parser.add_argument('--no_plotLT50', action='store_false', dest='plotLT50', help='add this argument to skip plotting the LT50 (3-2-1-0 scoring)')
    parser.add_argument('--representative', action='store', dest='rep', type=int, default=0, help='which number (specifying order/location) of the input files (1, 2, or 3, etc - based on indexing from 1) that is the representative to be used for I/LT50')
    parser.add_argument('--reportNum', action='store_false', dest='reportNum', help='add this argument to skip reporting the total number of worms in each concentration (summed across all input experiments), corresponding to Table 1 of the 2017 paper')
    parser.add_argument('--constrain1Hill', action='store', dest='hill1', type=float, default=-0.15, required=False, help='the constant/constrained Hill Slope for 1-0 scoring to be used when fewer than 3 replicates are provided')
    parser.add_argument('--constrain3Hill', action='store', dest='hill3', type=float, default=-1.5, required=False, help='the constant/constrained Hill Slope for 3-2-1-0 scoring to be used when fewer than 3 replicates are provided')
    parser.add_argument('--constrain2Hill', action='store', dest='hill2', type=float, default=-1.5, required=False, help='the constant/constrained Hill Slope for 2-1-0 scoring to be used when fewer than 3 replicates are provided')
    parser.add_argument('--constrain3Hill_bool', action='store_true', dest='constrain3Hill_bool', help='provide if and only if you want to constrain the hill slope for 3-2-1-0 scoring and then use the --constrain3Hill argument')
    parser.add_argument('--constrain2Hill_bool', action='store_true', dest='constrain2Hill_bool', help='provide if and only if you want to constrain the hill slope for 2-1-0 scoring and then use the --constrain2Hill argument')
    parser.add_argument('--constrain1Hill_bool', action='store_true', dest='constrain1Hill_bool', help='provide if and only if you want to constrain the hill slope for 1-0 scoring and then use the --constrain1Hill argument')
    parser.add_argument('--spline_k1', action='store', dest='spline_k1', type=int, default=3, required=False, help='the order of the first part of the spline smoothing if there is no _c50 fit')
    parser.add_argument('--spline_k2', action='store', dest='spline_k2', type=int, default=3, required=False, help='the order of the second part of the spline smoothing if there is no _c50 fit')
    parser.add_argument('--no_runNo3', action='store_false', dest='runNo3', help='add this argument to skip running additional analyses where the 3 & 2 scores are combined')
    parser.add_argument('--notDefaultC50', action='store_true', dest='notDefaultC50', help='provide this flag only if you want to provide and use all non-default initial parameters for C50 fits')
    parser.add_argument('--notDefTop3', action='store', dest='notDefTop3', type=float, default=100, help='not default Top initial parameter for IC50. Only provide if notDefaultC50 is also provided as True')
    parser.add_argument('--notDefTop1', action='store', dest='notDefTop1', type=float, default=100, help='not default Top initial parameter for LC50. Only provide if notDefaultC50 is also provided as True')
    parser.add_argument('--notDefBot3', action='store', dest='notDefBot3', type=float, default=0, help='not default Bottom initial parameter for IC50. Only provide if notDefaultC50 is also provided as True')
    parser.add_argument('--notDefBot1', action='store', dest='notDefBot1', type=float, default=0, help='not default Bottom initial parameter for LC50. Only provide if notDefaultC50 is also provided as True')
    parser.add_argument('--notDefic50', action='store', dest='notDefic50', type=float, default=1, help='not default ic50 initial parameter for IC50. Only provide if notDefaultC50 is also provided as True')
    parser.add_argument('--notDeflc50', action='store', dest='notDeflc50', type=float, default=10**1.5, help='not default lc50 initial parameter for LC50. Only provide if notDefaultC50 is also provided as True')
    parser.add_argument('--notDefHill3', action='store', dest='notDefHill3', type=float, default=-1, help='not default Hill Slope initial parameter for IC50. Only provide if notDefaultC50 is also provided as True')
    parser.add_argument('--notDefHill1', action='store', dest='notDefHill1', type=float, default=-1, help='not default Hill Slope initial parameter for LC50. Only provide if notDefaultC50 is also provided as True')
    parser.add_argument('--no_use_fitTop', action='store_false', dest='fitTopPlot_bool', help='use this flag if you want to plot with the 100/0 top/bottom instead of the fit top/bottom')

    parser.add_argument('--stats_compare_to_control', action='store_true', dest='stats_compare_to_control', help='use this flag if and only if you want to run chi square comparison to control; use --stats_concs and --stats_days to control what is compared and make sure to include at least one of --stats_mortality or --stats_inhibited')
    parser.add_argument('--stats_compare_to', action='store', nargs='+', dest='stats_compare_to', default=[0], help='all of the concentration indexes you want to be the expecation/compared to. 0 is 0 ug/mL. 1 is next highest concentration, etc.')
    parser.add_argument('--stats_days', action='store', nargs='+', type=int, default=[1,4,7], help='which days to compare the stats_concs to control')
    parser.add_argument('--stats_mortality', action='store_true', help='include if and only if you have included --stats_compare_to_control and want to test the difference in mortality')
    parser.add_argument('--stats_inhibited', action='store_true', help='include if and only if you have included --stats_compare_to_control and want to test the difference in number inhibited')
    #parser.add_argument('--random_seed', action='store', dest='random_seed', type=int, default=42)
    return parser

class WormAnalysis():
    def __init__(self, toAssess, expNames, drug, strain, stage, concUnits, molarmass, density, reportNum):
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

        self.consistent_num = True

        self.load_data(toAssess, expNames)
        self.find_mM(molarmass, density)
        if reportNum:
            self.find_numTotal()
        self.find_motility_index_score()
        self.find_mortality_score()

    def load_data(self, toAssess, expNames):
        self.num_experiments = len(toAssess)
        for i, file in enumerate(toAssess):
            df = pd.read_csv(file, sep='\t', header=0)
            if i == 0:
                '''find unique concentrations'''
                self.uniq_conc = np.unique(df.loc[:,'Concentration'])

                '''Make a data frame to store computed values in throughout the analysis; index is concentration'''
                index = []
                units = np.tile(self.concUnits_dict[self.concUnits], self.uniq_conc.shape[0])
                for conc, unit in zip(self.uniq_conc, units):
                    index.append(r'\textbf{%s %s}' % (conc,unit))
                self.df_tablec = pd.DataFrame(index=index)
                self.df_tablec.index.name='Concentration'

                '''set some dictionaries for mapping'''
                self.conc_index = {}
                self.index_to_conc = {}
                for k,c in enumerate(self.uniq_conc):
                    self.conc_index[float(c)] = k
                    self.index_to_conc[k] = float(c)

                '''find number of concentrations and days'''
                self.num_concentrations = len(self.uniq_conc)
                self.num_days = len(np.unique(df.loc[:,'Day']))

                '''Make a data frame to store computed values in throughout the analysis; index is Day'''
                days_arr = np.arange(self.num_days+1).astype(str)
                for k in range(days_arr.shape[0]):
                    days_arr[k] = r'\textbf{%s}' % days_arr[k]
                self.df_tabled = pd.DataFrame(index=days_arr)
                self.df_tabled.index.name='Day'
                self.df_tabled[r'\textbf{LC50}'] = np.tile('NC', days_arr.shape[0])
                self.df_tabled[r'\textbf{IC50}'] = np.tile('NC', days_arr.shape[0])

                '''initialize arrays for storing data'''
                self.scores3_by_well = np.zeros((self.num_concentrations*3, 4, self.num_days, self.num_experiments)) #num_concentrations*3 because concentrations for each experiment should be in triplicate, 4 because of 0-1-2-3 scoring, num_days, and num_experiments
                self.scores3_by_conc = np.zeros((self.num_concentrations, 4, self.num_days, self.num_experiments))

                '''set some more dictionaries for mapping; note that well index is base 0'''
                self.well_index_to_conc = {}
                self.conc_to_well_index = {}
                for well, conc in zip(np.array(df.loc[:, 'Well']).astype(np.int32)[0:self.num_concentrations*3], np.array(df.loc[:,'Concentration'])[0:self.num_concentrations*3]):
                    well_index = well-1
                    self.well_index_to_conc[well_index] = float(conc)
                    self.conc_to_well_index[float(conc)] = well_index

            '''fill arrays'''
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

                    to_add_val_0, to_add_val_1, to_add_val_2, to_add_val_3 = 0, 0, 0, 0
                    if len(fields[3]) > 0:
                        to_add_val_0 = int(fields[3])
                    if len(fields[4]) > 0:
                        to_add_val_1 = int(fields[4])
                    if len(fields[5]) > 0:
                        to_add_val_2 = int(fields[5])
                    if len(fields[6]) > 0:
                        to_add_val_3 = int(fields[6])
                    self.scores3_by_well[well_index, 0, day_index,i] = to_add_val_0 #store the num in well scoring 0
                    self.scores3_by_conc[conc_index, 0, day_index,i] += to_add_val_0
                    self.scores3_by_well[well_index, 1, day_index,i] = to_add_val_1 #store the num in well scoring 1
                    self.scores3_by_conc[conc_index, 1, day_index,i] += to_add_val_1
                    self.scores3_by_well[well_index, 2, day_index,i] = to_add_val_2 #store the num in well scoring 2
                    self.scores3_by_conc[conc_index, 2, day_index,i] += to_add_val_2
                    self.scores3_by_well[well_index, 3, day_index,i] = to_add_val_3 #store the num in well scoring 3
                    self.scores3_by_conc[conc_index, 3, day_index,i] += to_add_val_3

            num_trues = np.amax(np.unique(df.loc[:,'Well']).astype(np.int32))
            if num_trues != np.sum(df.loc[:,'numTotal_equal']):
                logging.warning('wells in {} do not have equal numbers of worms in a given well across the length of the experiment'.format(file))
                self.consistent_num = False
                fig, ax = plt.subplots(figsize=(10,10))

                max_total = 0
                for well_index in np.arange(self.scores3_by_well.shape[0]):
                    totals_by_day = np.sum(self.scores3_by_well[well_index, :, :, i], axis=0)
                    ax.plot(np.arange(self.num_days)+1, totals_by_day, marker='o', c=np.random.rand(3,), label="well {}".format(well_index + 1), clip_on=False)
                    if np.amax(totals_by_day) > max_total:
                        max_total = np.amax(totals_by_day)

                ax.set_ylim(0, max_total)
                ax.set_yticks(np.arange(0, max_total, 5))
                ax.set_yticklabels(np.arange(0, max_total, 5))
                ax.set_ylabel("Number of Worms")
                ax.set_xlabel("Day")
                ax.set_xticks(np.arange(self.num_days)+1)
                ax.set_xticklabels(np.arange(self.num_days)+1)
                ax.legend(bbox_to_anchor=[0.99, 0.5], loc='center left', ncol=2)
                plt.tight_layout()
                fig.savefig("well_num_worms_by_day_{}.pdf".format(expNames[i]))
                plt.close(fig)

    def reset_self_arrays(self, new_scores3_by_well, new_scores3_by_conc):
        '''Use this function with extreme caution'''
        '''Make sure that scores3_by_well is (num_concentrations*3, 4, num_days, num_experiments)
        and that scores3_by_conc is (num_concentrations, 4, num_days, num_experiments)'''
        logging.info("Following generation from the probabilistic model, resetting the 'input' data so we can repeat analysis with this generated data")
        self.old_scores3_by_well = np.copy(self.scores3_by_well)
        self.old_scores3_by_conc = np.copy(self.scores3_by_conc)
        self.scores3_by_well = new_scores3_by_well
        self.scores3_by_conc = new_scores3_by_conc
        self.find_motility_index_score()
        self.find_mortality_score()

    def find_mM(self, molarmass, density):
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
        mM_dict = {0: 1/molarmass/1e6*np.amax(self.uniq_conc)*1e3*1e3,
                   1: np.amax(self.uniq_conc)/1e3,
                   2: np.amax(self.uniq_conc)*1e3,
                   3: np.amax(self.uniq_conc),
                   4: np.amax(self.uniq_conc)/1e6,
                   5: 1/molarmass/1e9*np.amax(self.uniq_conc)*1e3*1e3,
                   6: 1/molarmass/1e3*np.amax(self.uniq_conc)*1e3*1e3,
                   7: 1/molarmass/1e3*np.amax(self.uniq_conc)*1e3*density*1e3*1e3,
                   8: 1/molarmass/1e3*np.amax(self.uniq_conc)*density*1e3*1e3}

        self.mM = '{0:.1f}'.format(mM_dict[self.concUnits])

    def find_numTotal(self):
        if self.consistent_num:
            total_nums = np.sum(self.scores3_by_conc[:,:,0,:].reshape(-1, 4*self.num_experiments), axis=1)
        else:
            logging.info('Note, since there were not a consistent number of worms in a given well for every day of the experiment, we find the total number of nematodes treated by the following steps'
                + '\nfirst, finding the day with the highest number of worms in each well for each experiment'
                + '\nsecond, summing the different experiments for each well'
                + '\nfinally, summing the triplicate wells for each concentration')
            total_nums = np.sum(np.sum(np.amax(np.sum(self.scores3_by_well, axis=1), axis=1), axis=1).reshape((self.num_concentrations, 3)), axis=1)
        self.total_num_in_well = np.amax(np.sum(self.scores3_by_well, axis=1), axis=1).astype(np.int32) #this is num_wells * num_experiments

        self.df_tablec[r'\textbf{Number of Nematodes Treated}'] = total_nums.astype(np.int32)
        logging.info("Added the total number of nematodes treated by concentration to the table which will be printed later. Column name is 'Number of Nematodes Treated'")


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
        right_side = ax.spines["right"]
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
        '''Specifies molarity of the highest concentration'''
        fig, ax = plt.subplots()

        days_arr = np.arange(self.num_days + 1)
        ax = self.format_plots(ax, title, xlabel, ylabel, ymin, ymax, ysep, days_arr)

        for i in range(toPlot.shape[0])[::-1]:
            if i == toPlot.shape[0]-1 and (self.concUnits in [0, 5, 6, 7, 8]):
                label = "{} {}\n= {} mM".format(self.uniq_conc[i], self.concUnits_dict[self.concUnits], self.mM)
            else:
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

    def driveLinePlots(self, plotLine3, plotLine1, isep, expNames, figname_base=''):
        '''Do we want to include the functionality to compute and plot stdev and sem? -- yes; do this'''
        if plotLine3:
            if isep:
                for i, exp in enumerate(expNames):
                    self.plotLineCurves(self.motility_index_scores_by_conc[:, :, i], 'isep_motility_{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, exp, figname_base), r"\textbf{%s %s on %s %s}" %(exp, self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$", "Motility Index Score", 0, 3, 1)
            reshaped_mid_by_well = self.motility_index_scores_by_well.reshape((self.num_concentrations, 3, self.num_days+1, self.num_experiments))
            motility_index_across_exp = np.zeros((self.num_concentrations, 3*self.num_experiments, self.num_days+1), dtype=np.float64)
            for j in range(self.num_experiments):
                motility_index_across_exp[:,j*3:(j*3)+3 ,:] = reshaped_mid_by_well[:,:,:,j]
            motility_avg_across_exp = np.mean(motility_index_across_exp, axis=1)
            self.plotLineCurves(motility_avg_across_exp, 'average_motility_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, figname_base), r"\textbf{%s on %s %s}" %(self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$", "Motility Index Score", 0, 3, 1)

        if plotLine1:
            if isep:
                for i, exp in enumerate(expNames):
                    self.plotLineCurves(self.mortality_scores_by_conc[:, :, i], 'isep_mortality_{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, exp, figname_base), r"\textbf{%s %s on %s %s}" %(exp, self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$" , "\% Alive", 0, 100, 25)
            reshaped_mort_by_well = self.mortality_scores_by_well.reshape((self.num_concentrations, 3, self.num_days+1, self.num_experiments))
            mortality_across_exp = np.zeros((self.num_concentrations, 3*self.num_experiments, self.num_days+1), dtype=np.float64)
            for j in range(self.num_experiments):
                mortality_across_exp[:,j*3:(j*3)+3, :] = reshaped_mort_by_well[:,:,:,j]
            mortality_avg_across_exp = np.mean(mortality_across_exp, axis=1)
            self.plotLineCurves(mortality_avg_across_exp, 'average_mortality_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, figname_base), r"\textbf{%s on %s %s}" %(self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$" , "\% Alive", 0, 100, 25)

    def plotSurvivalTime(self, inhibited, mortality, figname_base, plot_mortality=True, plot_motility=True, ysep=50, ymin=0, ymax=100, ylabel="\% Alive or \% Uninhibited", xlabel='Days'):
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
        if self.consistent_num:
            num_total = np.sum(toAnalyze_expSpec[:,:,0], axis=1).reshape((self.num_concentrations, 1))
        else:
            num_total = np.sum(toAnalyze_expSpec, axis=1)
        if motility:
            num_at_risk = toAnalyze_expSpec[:,3,:] #num at risk is only worms of score 3
            num_at_risk_corrected = np.minimum.accumulate(num_at_risk, axis=1)
            logging_value = r'\textbf{IT50}'
            logging_value2 = 'IT50'
        if mortality:
            num_at_risk = np.sum(toAnalyze_expSpec[:,1:,:], axis=1) #num at risk is any worm of score 1, 2, or 3
            num_at_risk_corrected = np.minimum.accumulate(num_at_risk, axis=1)
            logging_value = r'\textbf{LT50}'
            logging_value2 = 'LT50'
        toPopulate[:, 1:] = num_at_risk_corrected/num_total*100

        '''report day or U'''
        T50 = np.sum(toPopulate <= 50.0, axis=1)
        T50 = (self.num_days + 1 - T50).astype(str)
        bool_U = T50 == str(self.num_days + 1)
        T50[bool_U] = 'U'
        self.df_tablec[logging_value] = T50
        logging.info("Added the %s values to the table which will be printed later. Column name is '%s'" %(logging_value2, logging_value2))

        return toPopulate

    def driveSurvivalTimePlots(self, plotIT50, plotLT50, representative, expNames, figname_base=''):
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
            figname_base_spec = '{}_IT50_LT50_' + '{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, expNames[representative-1], figname_base)
            self.plotSurvivalTime(inhibited, mortality, figname_base_spec)
        else:
            if IT50:
                figname_base_spec = '{}_IT50' + '{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, expNames[representative-1], figname_base)
                self.plotSurvivalTime(inhibited, 0, figname_base_spec, plot_mortality=False)
            else:
                logging.warning('Why do you want to plot only the LT50? I suggest plotting them together. But here you go')
                figname_base_spec = '{}_LT50_' + '{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, expNames[representative-1], figname_base)
                self.plotSurvivalTime(0, mortality, figname_base_spec, plot_motility=False)
        logging.info('Completed Survival Analysis')

    def transform_X(self, X):
        '''set [0] to x0_val'''
        bool_0 = X == 0
        X[bool_0] += self.x0_val
        '''return log10 of concentrations'''
        return np.log10(X)

    def inhibitorResponse_equation(self, X, top, bottom, ic50, hillslope=-1):
        '''
        Y=Bottom + (Top-Bottom)/(1+10^((LogIC50-X)*HillSlope))
        where X is log10([]) with x=0 -> x=1e-6 (i.e. call self.transform_X(X))
        '''
        exponent = (np.log10(ic50)- self.transform_X(X))*hillslope
        equation = bottom + (top-bottom)/(1+(10**exponent))
        return equation

    def dose_response_sigmoid(self, X, X_mid, hill, bottom, top):
        '''
        For plotting found in def 114 here (https://gist.github.com/yannabraham/5f210fed773785d8b638)
        #Y = (c+(d-c)/(1+np.exp(b*(np.log(x)-np.log(e)))))
        '''
        #Y = (bottom+(top-bottom)/(1+np.exp(hill*(np.log(X)-np.log(X_mid)))))
        Y = (top + (bottom-top)/(1+np.exp(hill*(np.log(X) - np.log(X_mid)))))
        return Y

    def set_guess_params(self, avg1, avg3, hill1, hill3, default=True, num_exper='p3', not_default_3 = {}, not_default_1 = {}):
        if default:
            P0_30_top, P0_10_top = np.amax(avg3), np.amax(avg1)
            P0_30_bottom, P0_10_bottom = np.amin(avg3), np.amin(avg1)
            P0_30_ic50, P0_10_ic50 = 1, 10**1.5
            if num_exper == 'p3':
                P0_30_hill, P0_10_hill = -1, -1
                P0_30, P0_10 = [P0_30_top, P0_30_bottom, P0_30_ic50, P0_30_hill], [P0_10_top, P0_10_bottom, P0_10_ic50, P0_10_hill]
                return {'P0_30': P0_30,
                        'P0_10': P0_10,
                        'P0_30_top': P0_30_top,
                        'P0_30_bottom': P0_30_bottom,
                        'P0_30_ic50': P0_30_ic50,
                        'P0_30_hill': P0_30_hill,
                        'P0_10_top': P0_10_top,
                        'P0_10_bottom': P0_10_bottom,
                        'P0_10_ic50': P0_10_ic50,
                        'P0_10_hill': P0_10_hill}
            elif num_exper == 'l3':
                P0_30, P0_10 = [P0_30_top, P0_30_bottom, P0_30_ic50], [P0_10_top, P0_10_bottom, P0_10_ic50]
                return {'P0_30': P0_30,
                        'P0_10': P0_10,
                        'P0_30_top': P0_30_top,
                        'P0_30_bottom': P0_30_bottom,
                        'P0_30_ic50': P0_30_ic50,
                        'P0_30_hill': hill3,
                        'P0_10_top': P0_10_top,
                        'P0_10_bottom': P0_10_bottom,
                        'P0_10_ic50': P0_10_ic50,
                        'P0_10_hill': hill1}
        else:
            if num_exper == 'p3':
                P0_30 = [not_default_3['top'], not_default_3['bottom'], not_default_3['ic50'], not_default_3['hill']]
                P0_10 = [not_default_1['top'], not_default_1['bottom'], not_default_1['ic50'], not_default_1['hill']]
                return {'P0_30': P0_30,
                        'P0_10': P0_10,
                        'P0_30_top': not_default_3['top'],
                        'P0_30_bottom': not_default_3['bottom'],
                        'P0_30_ic50':  not_default_3['ic50'],
                        'P0_30_hill': not_default_3['hill'],
                        'P0_10_top': not_default_1['top'],
                        'P0_10_bottom': not_default_1['bottom'],
                        'P0_10_ic50': not_default_1['ic50'],
                        'P0_10_hill': not_default_1['hill']}
            elif num_exper == 'l3':
                P0_30 = [not_default_3['top'], not_default_3['bottom'], not_default_3['ic50']]
                P0_10 = [not_default_1['top'], not_default_1['bottom'], not_default_1['ic50']]
                return {'P0_30': P0_30,
                        'P0_10': P0_10,
                        'P0_30_top': not_default_3['top'],
                        'P0_30_bottom': not_default_3['bottom'],
                        'P0_30_ic50':  not_default_3['ic50'],
                        'P0_30_hill': hill3,
                        'P0_10_top': not_default_1['top'],
                        'P0_10_bottom': not_default_1['bottom'],
                        'P0_10_ic50': not_default_1['ic50'],
                        'P0_10_hill': hill1}

    def evaluate_no_fit(self, top, bottom, ic50, highest_conc):
        no_fit = False
        if top == bottom:
            no_fit = True
            logging.info('No fit possible since top and bottom are equal')
        elif bottom >= 50:
            no_fit = True
            logging.info('No fit possible since bottom is above 50 percent')
        elif ic50 > highest_conc:
            no_fit = True
            logging.info('No fit possible since reported ic50 is higher than highest_conc assayed')
        return no_fit

    def set_to_log_value(self, bottom, ic50, highest_conc, lowest_nz_conc, LC50=True):
        if LC50:
            logging_value = r'\textbf{LC50}'
            logging_value2 = 'LC50'
        else:
            logging_value = r'\textbf{IC50}'
            logging_value2 = 'IC50'
        logging_day = r'\textbf{%s}' % str(self.C_day)
        to_log = '{0:.3f}'.format(ic50)
        if bottom > 50 or ic50 > highest_conc:
            to_log = r'\textgreater' + '{}'.format(highest_conc)
        elif bottom < 50 and ic50 < lowest_nz_conc:
            to_log = r'\textless' + '{}'.format(lowest_nz_conc)
        self.df_tabled.loc[logging_day, logging_value] = to_log
        logging.info("Added the %s value to the table which will be printed later. Column name is '%s' and the day/row is '%s'" %(logging_value2, logging_value2, self.C_day))

    def plotIC(self, title, figname, concs_toAnalyze, averages, sems, curve_fit_ic50, curve_fit_hillslope = -1, curve_fit_top= 100, curve_fit_bottom=0, ylabel='\% Uninhibited', ysep=20, ymin=0, ymax=100):
        '''Let's try plotting with the Hill-Langmuir equation. Because plotting the inhibitor response curve isn't looking like the prism graphs'''
        conc_ticks = np.copy(self.uniq_conc)
        bool_0 = conc_ticks == 0
        conc_ticks[bool_0] = self.x0_val
        log_conc_ticks = np.log10(conc_ticks)

        concs = np.copy(concs_toAnalyze)
        bool_0 = concs == 0
        concs[bool_0] = self.x0_val
        log_concs = np.log10(concs)


        antilog_ax1_lim0 = self.x0_val
        antilog_ax1_lim1 = np.power(10, np.log10(self.x0_val)+1)
        antilog_ax2_lim0 = conc_ticks[1]
        antilog_ax2_lim1 = conc_ticks[-1]
        linspace_x1 = np.log10(np.linspace(antilog_ax1_lim0, antilog_ax1_lim1, 10e2))
        linspace_x2 = np.log10(np.linspace(antilog_ax2_lim0, antilog_ax2_lim1, 10e4))
        linspace_x1_antilog = np.power(np.tile(10, linspace_x1.shape[0]), linspace_x1)
        linspace_x2_antilog = np.power(np.tile(10, linspace_x2.shape[0]), linspace_x2)

        #curve1e = self.dose_response_sigmoid(linspace_x1_antilog, curve_fit_ic50, curve_fit_hillslope, curve_fit_top, curve_fit_bottom)
        curve1e = self.dose_response_sigmoid(linspace_x1_antilog, curve_fit_ic50, curve_fit_hillslope, curve_fit_bottom, curve_fit_top)
        #curve2e = self.dose_response_sigmoid(linspace_x2_antilog, curve_fit_ic50, curve_fit_hillslope, curve_fit_top, curve_fit_bottom)
        curve2e = self.dose_response_sigmoid(linspace_x2_antilog, curve_fit_ic50, curve_fit_hillslope, curve_fit_bottom, curve_fit_top)

        fig = plt.figure(constrained_layout=False)
        widths = [1, 8]
        gs = fig.add_gridspec(1, 2, width_ratios=widths, wspace=0.05)
        ax1 = fig.add_subplot(gs[0,0])
        ax2 = fig.add_subplot(gs[0,1])


        ax1.axhline(50, linestyle=':', color='black', clip_on=False)
        ax2.axhline(50, linestyle=':', color='black', clip_on=False)

        right_side1 = ax1.spines["right"]
        top_side1 = ax1.spines["top"]
        right_side2 = ax2.spines["right"]
        top_side2 = ax2.spines["top"]
        left_side2 = ax2.spines["left"]
        right_side1.set_visible(False)
        top_side1.set_visible(False)
        right_side2.set_visible(False)
        top_side2.set_visible(False)
        left_side2.set_visible(False)
        fig.suptitle(r'\textbf{%s}' %title)
        ax1.set_ylabel(r'\textbf{%s}' %ylabel)
        ax1.set_ylim(ymin, ymax, ysep)
        ax2.set_ylim(ymin, ymax, ysep)
        ax2.set_yticks([])
        ax2.set_yticklabels([])
        xlabel='Concentration (%s)' % self.concUnits_dict[self.concUnits]
        fig.text(0.5, 0.02, r'\textbf{%s}' %xlabel, ha='center', va='center')

        #mean values with SEM
        e = ax1.errorbar(log_concs[0], averages[0], yerr=sems[0], linestyle='None', marker='o', color='black', capsize=5, clip_on=False)
        for b in e[1]:
            b.set_clip_on(False)
        ax1.set_xlim(log_conc_ticks[0], log_conc_ticks[0]+1)
        ax1.set_xticks([log_conc_ticks[0], log_conc_ticks[0]+1])
        lead, power = str(self.x0_val).split("e-")
        ax1.set_xticklabels([r'$\mathrm{10^{-%s}}$' %str(int(power)), ' '])

        e = ax2.errorbar(log_concs[1:], averages[1:], yerr=sems[1:], linestyle='None', marker='o', color='black', capsize=5, clip_on=False)
        for b in e[1]:
            b.set_clip_on(False)
        ax2.set_xlim(log_conc_ticks[1], log_conc_ticks[-1])
        ax2.set_xticks(log_conc_ticks[1:])
        ticklabels = np.hstack((conc_ticks[1], conc_ticks[2:].astype(np.int32).astype(str)))
        ax2.set_xticklabels(ticklabels)

        ax1.plot(linspace_x1, curve1e, c='black', clip_on=False)

        ax2.plot(linspace_x2, curve2e, c='black', clip_on=False)

        ax1.plot((1,1), (1,1), color='black', clip_on=False) #bottom-left line
        ax2.plot((0,0), (1,1), color='black', clip_on=False) #bottom-right line

        fig.savefig(figname)
        plt.close(fig)
        logging.info('Plotted the figure {}'.format(figname))

    def plotIC_noFit(self, title, figname, concs_toAnalyze, averages, sems, spline_k1, spline_k2, ylabel='\% Uninhibited', ysep=20, ymin=0, ymax=100):
        '''Let's try plotting with the Hill-Langmuir equation. Because plotting the inhibitor response curve isn't looking like the prism graphs'''
        conc_ticks = np.copy(self.uniq_conc)
        bool_0 = conc_ticks == 0
        conc_ticks[bool_0] = self.x0_val
        log_conc_ticks = np.log10(conc_ticks)

        concs = np.copy(concs_toAnalyze)
        bool_0 = concs == 0
        concs[bool_0] = self.x0_val
        log_concs = np.log10(concs)

        fig = plt.figure(constrained_layout=False)
        widths = [1, 8]
        gs = fig.add_gridspec(1, 2, width_ratios=widths, wspace=0.05)
        ax1 = fig.add_subplot(gs[0,0])
        ax2 = fig.add_subplot(gs[0,1])


        ax1.axhline(50, linestyle=':', color='black', clip_on=False)
        ax2.axhline(50, linestyle=':', color='black', clip_on=False)

        right_side1 = ax1.spines["right"]
        top_side1 = ax1.spines["top"]
        right_side2 = ax2.spines["right"]
        top_side2 = ax2.spines["top"]
        left_side2 = ax2.spines["left"]
        right_side1.set_visible(False)
        top_side1.set_visible(False)
        right_side2.set_visible(False)
        top_side2.set_visible(False)
        left_side2.set_visible(False)
        fig.suptitle(r'\textbf{%s}' %title)
        ax1.set_ylabel(r'\textbf{%s}' %ylabel)
        ax1.set_ylim(ymin, ymax, ysep)
        ax2.set_ylim(ymin, ymax, ysep)
        ax2.set_yticks([])
        ax2.set_yticklabels([])
        xlabel='Concentration (%s)' % self.concUnits_dict[self.concUnits]
        fig.text(0.5, 0.02, r'\textbf{%s}' %xlabel, ha='center', va='center')

        #mean values with SEM
        e = ax1.errorbar(log_concs[0], averages[0], yerr=sems[0], marker='o', fmt='none', linestyle=None, color='black', capsize=5, clip_on=False)
        for b in e[1]:
            b.set_clip_on(False)

        x_s_1 = np.linspace(log_conc_ticks[0], log_conc_ticks[-1], 450)
        spl_1 = make_interp_spline(log_concs, averages, k=spline_k1)
        power_smooth_1 = spl_1(x_s_1)
        s = ax1.scatter(log_conc_ticks[0], averages[0], c='black')
        ax1.plot(x_s_1, power_smooth_1, c='black')

        s.set_clip_on(False)

        ax1.set_xlim(log_conc_ticks[0], log_conc_ticks[0]+1)
        ax1.set_xticks([log_conc_ticks[0], log_conc_ticks[0]+1])
        lead, power = str(self.x0_val).split("e-")
        ax1.set_xticklabels([r'$\mathrm{10^{-%s}}$' %str(int(power)), ' '])

        e = ax2.errorbar(log_concs[1:], averages[1:], yerr=sems[1:], marker='o', fmt='none', linestyle=None, color='black', capsize=5, clip_on=False)
        for b in e[1]:
            b.set_clip_on(False)

        x_s_2 = np.linspace(log_concs[1], log_conc_ticks[-1], 300)
        spl_2 = make_interp_spline(log_concs[1:], averages[1:], k=spline_k2)
        power_smooth_2 = spl_2(x_s_2)
        idx = np.argwhere(np.diff(np.sign(power_smooth_2 - 50))).flatten()
        print(np.power(10, x_s_2[idx]), flush=True)
        s = ax2.scatter(log_conc_ticks[1:], averages[1:], c='black')
        ax2.plot(x_s_2, power_smooth_2, c='black', clip_on=False)

        s.set_clip_on(False)

        ax2.set_xlim(log_conc_ticks[1], log_conc_ticks[-1])
        ax2.set_xticks(log_conc_ticks[1:])
        ticklabels = np.hstack((conc_ticks[1], conc_ticks[2:].astype(np.int32).astype(str)))
        ax2.set_xticklabels(ticklabels)

        ax1.plot((1,1), (1,1), color='black', clip_on=False) #bottom-left line
        ax2.plot((0,0), (1,1), color='black', clip_on=False) #bottom-right line

        fig.savefig(figname)
        plt.close(fig)
        logging.info('Plotted the figure {}'.format(figname))

    def driveIC(self, plotIC50, plotLC50, C_day, x0_val, hill1, hill3, spline_k1, spline_k2, fitTopPlot_bool, hill3bool, hill1bool, default=True, not_default_3 = {}, not_default_1 = {}, figname_base=''):
        #Look at each well self.scores3_by_well
        #recall self.scores3_by_well = np.zeros((self.num_concentrations*3, 4, self.num_days, self.num_experiments))
        #Use these to go from well to conc etc
        #self.well_index_to_conc = {}
        #self.conc_to_well_index = {}
        self.C_day = C_day
        self.x0_val = x0_val
        totals = np.sum(self.scores3_by_well[:,:,self.C_day-1,:], axis=1)

        uninhibited3 = self.scores3_by_well[:,3,self.C_day-1,:]/totals*100
        uninhibited1 = np.sum(self.scores3_by_well[:,1:,self.C_day-1,:], axis=1)/totals*100

        '''find averages and SEMs (note: we plot SEMs because we want to look at the fit of the data to the curve fit model)'''
        toAvg_3 = np.zeros((self.num_concentrations, 3*self.num_experiments))
        toAvg_1 = np.zeros((self.num_concentrations, 3*self.num_experiments))
        for i in range(uninhibited3.shape[0]):
            toAvg_3[self.conc_index[self.well_index_to_conc[i]], ((i%3) * self.num_experiments) : (((i%3)+1) * self.num_experiments)] = uninhibited3[i].copy()
            toAvg_1[self.conc_index[self.well_index_to_conc[i]], ((i%3) * self.num_experiments) : (((i%3)+1) * self.num_experiments)] = uninhibited1[i].copy()

        sem3 = sem(toAvg_3, axis=1)
        avg3 = np.average(toAvg_3, axis=1)

        sem1 = sem(toAvg_1, axis=1)
        avg1 = np.average(toAvg_1, axis=1)

        conc_X = np.tile(np.array([self.well_index_to_conc[x] for x in np.arange(self.num_concentrations*3)]).reshape(-1, 1), (1, self.num_experiments))

        if self.num_experiments >= 3 and default and not hill3bool:
            param_dict = self.set_guess_params(avg1, avg3, hill1, hill3)

        elif self.num_experiments >= 3 and default and hill3bool:
            param_dict = self.set_guess_params(avg1, avg3, hill1, hill3, num_exper='l3')

        elif self.num_experiments < 3 and default and not hill3bool:
            param_dict = self.set_guess_params(avg1, avg3, hill1, hill3, num_exper='l3')

        elif self.num_experiments >= 3 and not default:
            param_dict = self.set_guess_params(avg1, avg3, hill1, hill3, default=False, not_default_3=not_default_3, not_default_1=not_default_1 )

        elif self.num_experiments < 3 and not default:
            param_dict = self.set_guess_params(avg1, avg3, hill1, hill3, default=False, num_exper='l3', not_default_3=not_default_3, not_default_1=not_default_1)

        P0_30 = param_dict['P0_30']
        P0_10 = param_dict['P0_10']
        P0_30_top = param_dict['P0_30_top']
        P0_30_bottom = param_dict['P0_30_bottom']
        P0_30_ic50 = param_dict['P0_30_ic50']
        P0_30_hill = param_dict['P0_30_hill']
        P0_10_top = param_dict['P0_10_top']
        P0_10_bottom = param_dict['P0_10_bottom']
        P0_10_ic50 = param_dict['P0_10_ic50']
        P0_10_hill = param_dict['P0_10_hill']

        lowest_nz_conc = np.sort(self.uniq_conc)[1]
        highest_conc = np.amax(self.uniq_conc)

        if self.num_experiments >=3 and not hill3bool:
            if plotLC50:
                logging.info('Running Levenberg-Marquardt Algorithm Scipy Curve Fitting for 1-0 scoring using a max number of function evaluations of {}. Initial values are the following.\nINITIAL Top:\t{}\nINITIAL Bottom:\t{}\nINITIAL IC50:\t{}\nINITIAL HillSlope:\t{}'.format(int(1e6), P0_10_top, P0_10_bottom, P0_10_ic50, P0_10_hill))
                popt, popc = curve_fit(self.inhibitorResponse_equation, conc_X.flatten(), uninhibited1.flatten(), p0=P0_10, method='lm', maxfev=int(1e6))
                top_10, bottom_10, ic50_10, hill_10 = popt[0], popt[1], popt[2], popt[3]
                logging.info('Returned lm fit for 1-0 scoring.\nFIT Top:\t{}\nFIT Bottom:\t{}\nFIT IC50:\t{}\nFIT HillSlope:\t{}'.format(top_10, bottom_10, ic50_10, hill_10))
                no_fit_bool = self.evaluate_no_fit(top_10, bottom_10, ic50_10, highest_conc)
                self.set_to_log_value(bottom_10, ic50_10, highest_conc, lowest_nz_conc)
                if no_fit_bool:
                    self.plotIC_noFit(r'$\mathrm{LC_{50}}$'+ ' {} on {} {} Day {}'.format(self.drug, self.stage, self.strain, self.C_day), 'LC50_{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, self.C_day, figname_base), self.uniq_conc, avg1, sem1, spline_k1, spline_k2)
                else:
                    if fitTopPlot_bool:
                        self.plotIC(r'$\mathrm{LC_{50}}$' + ' {} on {} {} Day {}'.format(self.drug, self.stage, self.strain, self.C_day), 'LC50_{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, self.C_day, figname_base), self.uniq_conc, avg1, sem1, ic50_10, curve_fit_hillslope=hill_10, curve_fit_top=top_10, curve_fit_bottom = bottom_10)
                    else:
                        self.plotIC(r'$\mathrm{LC_{50}}$' + ' {} on {} {} Day {}'.format(self.drug, self.stage, self.strain, self.C_day), 'LC50_{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, self.C_day, figname_base), self.uniq_conc, avg1, sem1, ic50_10, cruve_fit_hillslope=hill_10)
            if plotIC50:
                logging.info('Running Levenberg-Marquardt Algorithm Scipy Curve Fitting for 3-2-1-0 scoring using the default max number of function evaluations. Initial values are the following.\nINITIAL Top:\t{}\nINITIAL Bottom:\t{}\nINITIAL IC50:\t{}\nINITIAL HillSlope:\t{}'.format(P0_30_top, P0_30_bottom, P0_30_ic50, P0_30_hill))
                popt2, popc2 = curve_fit(self.inhibitorResponse_equation, conc_X.flatten(), uninhibited3.flatten(), p0=P0_30, method='lm', maxfev=int(1e6))
                top_30, bottom_30, ic50_30, hill_30 = popt2[0], popt2[1], popt2[2], popt2[3]
                logging.info('Returned lm fit for 3-2-1-0 scoring.\nFIT Top:\t{}\nFIT Bottom:\t{}\nFIT IC50:\t{}\nFIT HillSlope:\t{}'.format(top_30, bottom_30, ic50_30, hill_30))
                no_fit_bool = self.evaluate_no_fit(top_30, bottom_30, ic50_30, highest_conc)
                self.set_to_log_value(bottom_30, ic50_30, highest_conc, lowest_nz_conc, LC50=False)
                if no_fit_bool:
                    self.plotIC_noFit(r'$\mathrm{IC_{50}}$' + ' {} on {} {} Day {}'.format(self.drug, self.stage, self.strain, self.C_day), 'IC50_{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, self.C_day, figname_base), self.uniq_conc, avg3, sem3, spline_k1, spline_k2)
                else:
                    if fitTopPlot_bool:
                        self.plotIC(r'$\mathrm{IC_{50}}$' + ' {} on {} {} Day {}'.format(self.drug, self.stage, self.strain, self.C_day), 'IC50_{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, self.C_day, figname_base), self.uniq_conc, avg3, sem3, ic50_30, curve_fit_hillslope=hill_30, curve_fit_top=top_30, curve_fit_bottom=bottom_30)
                    else:
                        self.plotIC(r'$\mathrm{IC_{50}}$' + ' {} on {} {} Day {}'.format(self.drug, self.stage, self.strain, self.C_day), 'IC50_{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, self.C_day, figname_base), self.uniq_conc, avg3, sem3, ic50_30, curve_fit_hillslope=hill_30)
            logging.info('Completed Non-linear Regression for Inhibition Response Analysis')

        else:
            if plotLC50:
                logging.info('Running Levenberg-Marquardt Algorithm Scipy Curve Fitting for 1-0 scoring using a max number of function evaluations of {} and a constant hill slope of {}. Initial values are the following.\nINITIAL Top:\t{}\nINITIAL Bottom:\t{}\nINITIAL IC50:\t{}'.format(int(1e6), P0_10_hill, P0_10_top, P0_10_bottom, P0_10_ic50))
                popt, popc = curve_fit(self.inhibitorResponse_equation, conc_X.flatten(), uninhibited1.flatten(), p0=P0_10, method='lm', maxfev=int(1e6))
                top_10, bottom_10, ic50_10 = popt[0], popt[1], popt[2]
                logging.info('Returned lm fit for 1-0 scoring.\nFIT Top:\t{}\nFIT Bottom:\t{}\nFIT IC50:\t{}'.format(top_10, bottom_10, ic50_10))
                no_fit_bool = self.evaluate_no_fit(top_10, bottom_10, ic50_10, highest_conc)
                self.set_to_log_value(bottom_10, ic50_10, highest_conc, lowest_nz_conc)
                if no_fit_bool:
                    self.plotIC_noFit(r'$\mathrm{LC_{50}}$' + ' {} on {} {} Day {}'.format(self.drug, self.stage, self.strain, self.C_day), 'LC50_{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, self.C_day, figname_base), self.uniq_conc, avg1, sem1, spline_k1, spline_k2)
                else:
                    if fitTopPlot_bool:
                        self.plotIC(r'$\mathrm{LC_{50}}$' + ' {} on {} {} Day {}'.format(self.drug, self.stage, self.strain, self.C_day), 'LC50_{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, self.C_day, figname_base), self.uniq_conc, avg1, sem1, ic50_10, curve_fit_hillslope=P0_10_hill, curve_fit_top=top_10, curve_fit_bottom=bottom_10)
                    else:
                        self.plotIC(r'$\mathrm{LC_{50}}$' + ' {} on {} {} Day {}'.format(self.drug, self.stage, self.strain, self.C_day), 'LC50_{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, self.C_day, figname_base), self.uniq_conc, avg1, sem1, ic50_10, curve_fit_hillslope=P0_10_hill)
            if plotIC50:
                logging.info('Running Levenberg-Marquardt Algorithm Scipy Curve Fitting for 3-2-1-0 scoring using the default max number of function evaluations and a constant hill slope of {}. Initial values are the following.\nINITIAL Top:\t{}\nINITIAL Bottom:\t{}\nINITIAL IC50:\t{}'.format(P0_30_hill, P0_30_top, P0_30_bottom, P0_30_ic50))
                popt2, popc2 = curve_fit(self.inhibitorResponse_equation, conc_X.flatten(), uninhibited3.flatten(), p0=P0_30, method='lm', maxfev=int(1e6))
                top_30, bottom_30, ic50_30 = popt2[0], popt2[1], popt2[2]
                logging.info('Returned lm fit for 3-2-1-0 scoring.\nFIT Top:\t{}\nFIT Bottom:\t{}\nFIT IC50:\t{}'.format(top_30, bottom_30, ic50_30))
                no_fit_bool = self.evaluate_no_fit(top_30, bottom_30, ic50_30, highest_conc)
                self.set_to_log_value(bottom_30, ic50_30, highest_conc, lowest_nz_conc, LC50=False)
                if no_fit_bool:
                    self.plotIC_noFit(r'$\mathrm{IC_{50}}$' + ' {} on {} {} Day {}'.format(self.drug, self.stage, self.strain, self.C_day), 'IC50_{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, self.C_day, figname_base), self.uniq_conc, avg3, sem3, spline_k1, spline_k2)
                else:
                    if fitTopPlot_bool:
                        self.plotIC(r'$\mathrm{IC_{50}}$' + ' {} on {} {} Day {}'.format(self.drug, self.stage, self.strain, self.C_day), 'IC50_{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, self.C_day, figname_base), self.uniq_conc, avg3, sem3, ic50_30, curve_fit_hillslope=P0_30_hill, curve_fit_top=top_30, curve_fit_bottom=bottom_30)
                    else:
                        self.plotIC(r'$\mathrm{IC_{50}}$' + ' {} on {} {} Day {}'.format(self.drug, self.stage, self.strain, self.C_day), 'IC50_{}_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, self.C_day, figname_base), self.uniq_conc, avg3, sem3, ic50_30, curve_fit_hillslope=P0_30_hill)
            logging.info('Completed Non-linear Regression for Inhibition Response Analysis but with constrained hill slope therefore the fit is likely less than ideal')

    def reportTable(self, rep_exp, reportNum, plotIT50, plotLT50, plotIC50, plotLC50, figname_base=''):
        '''write pandas dataframes of computed values [(number of worms treated, IT50, LT50) and (IC50, LC50) to a pdf file with table(s)'''
        if (reportNum or plotIT50 or plotLT50) and (plotIC50 or plotLC50):
            fig, (ax1, ax2) = plt.subplots(2,1)
            #hide axes
            fig.patch.set_visible(False)
            ax1.axis('off')
            ax2.axis('off')
            ax1.axis('tight')
            ax2.axis('tight')
            #write table
            ax1.table(cellText=self.df_tablec.values, colLabels=self.df_tablec.columns, rowLabels=self.df_tablec.index, cellLoc = 'center', loc='center')
            ax2.table(cellText=self.df_tabled.values, colLabels=self.df_tabled.columns, rowLabels=self.df_tabled.index, cellLoc = 'center', loc='center')
            filename='table_{}_{}_{}_t50repexp_{}_c50day_{}{}.pdf'.format(self.drug, self.stage, self.strain, rep_exp, self.C_day, figname_base)

        elif (plotIC50 or plotLC50) and not (reportNum or plotIT50 or plotLT50):
            fig, ax = plt.subplots()
            fig.patch.set_visible(False)
            ax.axis('off')
            ax.axis('tight')
            ax.table(cellText=self.df_tabled.values, colLabels=self.df_tabled.columns, rowLabels=self.df_tabled.index, cellLoc = 'center', loc='center')
            filename='table_{}_{}_{}_c50day_{}{}.pdf'.format(self.drug, self.stage, self.strain, self.C_day, figname_base)

        elif (reportNum) and not (plotIC50 or plotLC50 or plotIT50 or plotLT50):
            fig, ax = plt.subplots()
            fig.patch.set_visible(False)
            ax.axis('off')
            ax.axis('tight')
            ax.table(cellText=self.df_tablec.values, colLabels=self.df_tablec.columns, rowLabels=self.df_tablec.index, cellLoc = 'center', loc='center')
            filename = 'table_{}_{}_{}{}.pdf'.format(self.drug, self.stage, self.strain, figname_base)

        elif (plotLT50 or plotIT50) and not (plotIC50 or plotLC50):
            fig, ax = plt.subplots()
            fig.patch.set_visible(False)
            ax.axis('off')
            ax.axis('tight')
            ax.table(cellText=self.df_tablec.values, colLabels=self.df_tablec.columns, rowLabels=self.df_tablec.index, cellLoc = 'center', loc='center')
            filename = 'table_{}_{}_{}_t50repexp_{}{}.pdf'.format(self.drug, self.stage, self.strain, rep_exp, figname_base)

        fig.tight_layout()
        fig.savefig(filename, format='pdf')
        plt.close(fig)
        logging.info('Wrote table(s) of computed values to {}'.format(filename))

    def find_num_mortality(self, vec_of_scores):
        return (vec_of_scores[0])

    def find_opp_nm(self, vec_of_scores):
        return (np.sum(vec_of_scores[1:]))

    def find_num_inhibited(self, vec_of_scores):
        return (np.sum(vec_of_scores[0:3]))

    def find_opp_ni(self, vec_of_scores):
        return(vec_of_scores[-1])

    def make_mor_arr(self, vec_of_scores):
        return (np.hstack((self.find_num_mortality(vec_of_scores), self.find_opp_nm(vec_of_scores))))

    def make_inh_arr(self, vec_of_scores):
        return (np.hstack((self.find_num_inhibited(vec_of_scores), self.find_opp_ni(vec_of_scores))))

    def driveStats(self, compare_to_conc, compare_days, mortality_bool, motility_bool):
        for compare_to in compare_to_conc:
            compare_to_index = self.conc_index[compare_to]
            compare_with_concs = self.uniq_conc[self.uniq_conc != compare_to]
            if motility_bool:
                df_table_i = pd.DataFrame(index=compare_with_concs)
                df_table_i.index.name = "Concentration (ug/mL)"
                df_table_i_bools = pd.DataFrame(index=compare_with_concs)
                df_table_i_bools.index.name = "Concentration (ug/mL)"
            if mortality_bool:
                df_table_m = pd.DataFrame(index=compare_with_concs)
                df_table_m.index.name = "Concentration (ug/mL)"
                df_table_m_bools = pd.DataFrame(index=compare_with_concs)
                df_table_m_bools.index.name = "Concentration (ug/mL)"
            for compare_day in compare_days:
                if motility_bool:
                    compare_to_inh = self.make_inh_arr(np.sum(self.scores3_by_conc[compare_to_index, :, compare_day -1, :], axis=1))
                    pvals_i = np.array([])
                    pvals_i_bools = np.array([])
                if mortality_bool:
                    compare_to_mor = self.make_mor_arr(np.sum(self.scores3_by_conc[compare_to_index, :, compare_day -1, :], axis=1))
                    pvals_m = np.array([])
                    pvals_m_bools = np.array([])
                for compare_with in compare_with_concs:
                    compare_with_index = self.conc_index[compare_with]
                    if motility_bool:
                        compare_with_inh = self.make_inh_arr(np.sum(self.scores3_by_conc[compare_with_index, :, compare_day -1, :], axis=1))
                        if np.sum(compare_with_inh > 5) == 2 and np.sum(compare_to_inh > 5) == 2:
                            #print(compare_with_inh, compare_to_inh)
                            chisq, pval_i = chisquare(compare_with_inh, f_exp = compare_to_inh)
                            pval_bool = pval_i < 0.05
                            pval_i = format(pval_i, ".2e")

                        else:
                            pval_i = "U"
                            pval_bool = "NC"
                        pvals_i = np.hstack((pvals_i, pval_i))
                        pvals_i_bools = np.hstack((pvals_i_bools, pval_bool))
                    if mortality_bool:
                        compare_with_mor = self.make_mor_arr(np.sum(self.scores3_by_conc[compare_with_index, :, compare_day -1, :], axis=1))
                        if np.sum(compare_with_mor > 5) == 2 and np.sum(compare_to_mor > 5) == 2:
                            #print(compare_with_mor, compare_to_mor)
                            chisq, pval_m = chisquare(compare_with_mor, f_exp = compare_to_mor)
                            pval_bool = pval_m < 0.05
                            pval_m = format(pval_m, ".2e")
                        else:
                            pval_m = "U"
                            pval_bool = "NC"
                        pvals_m = np.hstack((pvals_m, pval_m))
                        pvals_m_bools = np.hstack((pvals_m_bools, pval_bool))


                if motility_bool:
                    df_table_i["Day {}".format(compare_day)] = pvals_i
                    df_table_i_bools["Day {}".format(compare_day)] = pvals_i_bools
                    df_table_i_bools["Day {}".format(compare_day)][df_table_i_bools["Day {}".format(compare_day)] != "NC"] = df_table_i_bools["Day {}".format(compare_day)][df_table_i_bools["Day {}".format(compare_day)] != "NC"].astype('bool')
                if mortality_bool:
                    df_table_m["Day {}".format(compare_day)] = pvals_m
                    df_table_m_bools["Day {}".format(compare_day)] = pvals_m_bools
                    df_table_m_bools["Day {}".format(compare_day)][df_table_m_bools["Day {}".format(compare_day)] != "NC"] = df_table_m_bools["Day {}".format(compare_day)][df_table_m_bools["Day {}".format(compare_day)] != "NC"].astype('bool')
            if motility_bool:
                fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
                fig.suptitle("chi-square analysis of fraction inhibited with expectation from {} ug/mL".format(compare_to))
                fig.patch.set_visible(False)
                for ax in [ax1, ax2]:
                    ax.axis('off')
                    ax.axis('tight')
                ax1.set_title("p-value")
                ax1.table(cellText=df_table_i.values, colLabels=df_table_i.columns, rowLabels = df_table_i.index, cellLoc='center', loc='center')
                ax2.set_title('p-value' + r'$< 0.05$')
                ax2.table(cellText=df_table_i_bools.values, colLabels=df_table_i.columns, rowLabels = df_table_i.index, cellLoc='center', loc='center')
                fig.tight_layout()
                fig.savefig("pval_table_inhibition_expectation{}_{}_{}_{}.pdf".format(compare_to, self.drug, self.stage, self.strain))
                plt.close(fig)
            if mortality_bool:
                fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
                fig.suptitle("chi-square analysis of fraction mortality with expectation from {} ug/mL".format(compare_to))
                fig.patch.set_visible(False)
                for ax in [ax1, ax2]:
                    ax.axis('off')
                    ax.axis('tight')
                ax.set_title("p-value")
                ax1.table(cellText=df_table_m.values, colLabels=df_table_m.columns, rowLabels = df_table_m.index, cellLoc='center', loc='center')
                ax2.set_title('p-value' + r'$ < 0.05$')
                ax2.table(cellText=df_table_m_bools.values, colLabels=df_table_i.columns, rowLabels = df_table_i.index, cellLoc='center', loc='center')
                fig.tight_layout()
                fig.savefig("pval_table_mortality_expectation{}_{}_{}_{}.pdf".format(compare_to, self.drug, self.stage, self.strain))
                plt.close(fig)


main()
