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

def main():
    parser = generate_parser()
    args = parser.parse_args()
    analysis_instance = WormAnalysis(args.toAssess, args.logfile)
    if args.plotIC50 or args.plotLC50:
        analysis_instance.driveIC(args.plotIC50, args.plotLC50, args.C_day, args.x0_val)


def generate_parser():
    parser = ap.ArgumentParser(description='C elegans analysis')
    parser.add_argument('--toAssess', action='store', dest='toAssess', nargs='+', type=str, required=True, help='the exported files to assess')
    parser.add_argument('--plotLine3', action='store', dest='plotLine3', type=bool, default=True, help='whether to plot the daily motility response by concentration')
    parser.add_argument('--plotLine1', action='store', dest='plotLine1', type=bool, default=True, help='whether to plot the daily lethality response by concentration')
    parser.add_argument('--plotIC50', action='store', dest='plotIC50', type=bool, default=True, help='whether to plot the IC50 (3-2-1-0 scoring)')
    parser.add_argument('--plotLC50', action='store', dest='plotLC50', type=bool, default=True, help='whether to plot the LC50 (1-0 scoring)')
    parser.add_argument('--C_day', action='store', dest='C_day', type=int, default=4, help='the day (index from 1) to assess for inhibitory and/or lethal concentration')
    parser.add_argument('--x0_value', action='store', dest='x0_val', type=float, default=1e-6, help='value to replace the x=0 [] with when transforming x')
    parser.add_argument('--plotIT50', action='store', dest='plotIT50', type=bool, default=True, help='whether to plot the IT50 (3-2-1-0 scoring)')
    parser.add_argument('--plotLT50', action='store', dest='plotLT50', type=bool, default=True, help='whether to plot the LT50 (3-2-1-0 scoring)')
    parser.add_argument('--representative', action='store', dest='rep', type=int, default=0, help='which number (specifying order/location) of the input files (1, 2, or 3, etc) that is the representative to be used for I/LT50')
    parser.add_argument('--logfile', action='store', dest='logfile', type=str, default='worm_analysis_{}.txt'.format(datetime.datetime.now()), help='the logfile to store information and any warning messages')
    return parser

class WormAnalysis():
    def __init__(self, toAssess, logfile):
        logging.basicConfig(filename=logfile, level=logging.INFO)
        self.load_data(toAssess)
        self.find_motility_index_score()
        self.find_mortality_score()

    def load_data(self, toAssess):
        num_experiments = len(toAssess)
        for i, file in enumerate(toAssess):
            df = pd.read_csv(file, sep='\t', header=0)
            if i == 0:
                uniq_conc = np.unique(df.loc[:,'Concentration'])
                self.conc_index = {}
                self.index_to_conc = {}
                for k,c in enumerate(uniq_conc):
                    self.conc_index[float(c)] = k
                    self.index_to_conc[k] = float(c)
                num_concentrations = len(uniq_conc)
                num_days = len(np.unique(df.loc[:,'Day']))
                self.scores3_by_well = np.zeros((num_concentrations*3, 4, num_days, num_experiments)) #num_concentrations*3 because concentrations for each experiment should be in triplicate, 4 because of 0-1-2-3 scoring, num_days, and num_experiments
                self.scores3_by_conc = np.zeros((num_concentrations, 4, num_days, num_experiments))
                self.well_index_to_conc = {}
                self.conc_to_well_index = {}
                for well, conc in zip(np.array(df.loc[:, 'Well']).astype(np.int32)[0:num_concentrations*3], np.array(df.loc[:,'Concentration'])[0:num_concentrations*3]):
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

    def find_motility_index_score(self):
        return 0

    def find_mortality_score(self):
        return 0

    def transform_X(self, X):
        return 0

    def inhibitorResponse_equation(self, X, top, bottom, ic50, hillslope):
        exponent = (np.log10(ic50)- self.transform_X(X))*hillslope
        equation = bottom + (top-bottom)/(1+(10**exponent))
        return equation

    def driveIC(self, plotIC50, plotLC50, C_day, x0_val):
        self.C_day = C_day
        self.x0_val = x0_val

main()
