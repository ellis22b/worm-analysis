#!/usr/bin/env python3

import numpy as np
import logging
import datetime
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt



class WormAnalysis_probModel():
    def __init__(self, drug, strain, stage, num_days, num_experiments, num_concentrations, scores3_by_well, scores3_by_conc, total_nums, total_num_in_well, conc_index, well_index_to_conc, random_seed):
        self.drug = drug
        self.strain = strain
        self.stage = stage
        self.num_days = num_days
        self.num_experiments = num_experiments
        self.num_concentrations = num_concentrations
        self.scores3_by_well = scores3_by_well
        self.scores3_by_conc = scores3_by_conc
        self.total_nums = total_nums
        self.total_num_in_well = total_num_in_well
        self.conc_index = conc_index
        self.well_index_to_conc = well_index_to_conc
        self.random_seed = random_seed

        logging.basicConfig(filename=new_log_file, level=logging.INFO, filemode='w', format='%(name)s - %(levelname)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')

    def run(self):
        params = self.learn_PWM_params() #this is number of wells, number of scores, number of days + 1, number of experiments
        self.plotInformationContent(params)
        sequence_paths = self.generate_sequences(params, self.total_num_in_well)
        new_scores3_by_well, new_scores3_by_conc = self.recreate_scores3(sequence_paths)
        return(new_scores3_by_well, new_scores3_by_conc)

    def get_num_treated(self):
        return np.array(self.df_tablec.loc[:,r'\textbf{Number of Nematodes Treated}'])

    def learn_PWM_params(self):
        '''Here we'll use the number that are actually observed at each day to learn the params;
        when we generate based on these parameters later, we'll generate the minimum number that should have been observed each day'''
        logging.info('Because of incomplete data, learning a probabilistic model to rerun all of the analysis with generated sequences')
        logging.info('learning PWM parameters')
        params = np.zeros((self.num_concentrations *3, 4, self.num_days+1, self.num_experiments))

        '''use self.scores3_by_well[well_i, score, day_index, exp_num] to find numerator and denominator'''
        for j in range(self.num_experiments):
            params[:,:,0, j] = np.tile([0, 0, 0, 1], (self.num_concentrations*3, 1))#on day 0, only scores of 3
            totals = np.sum(self.scores3_by_well[:,:,:,j], axis=1).reshape((self.num_concentrations*3, 1, self.num_days))
            params[:, :, 1:, j] = self.scores3_by_well[:,:,:,j]/totals
        return params

    def plotInformationContent(self, PWM):
        '''The information content (y-axis) of position i is given by
            R_i = log_2 (4) - (H_i + e_n)
            where H_i is the uncertainty of position i given by
            H_i = - sum_{b=0 to 3} f_{b,i} * log_2 (f_{b,i})
            where f_{b,i} is the relative frequency of score b at position i
            and e_n is the small-sample correction for an alignment of n scores
            e_n = 1 / ln (2) * (s-1) / 2n
            where s = 4 and n is the number of sequences in the alignment
            The height of score a in column i is given by
            height = f_{b,i} * R_i'''

        fig, axes = plt.subplots(nrows = np.ceil(num_concentrations * 3 /2), ncols = 6)
        xlabel = "day"
        ylabel = "bits"
        for i, ax in enumerate(axes):
            right_side = ax.spines["right"]
            top_side = ax.spines["top"]
            right_side.set_visible(False)
            top_side.set_visible(False)
            ax.set_title("Well {}".format(i), fontsize=12)
            ax.set_xlabel(r'\textbf{%s}' %xlabel, fontsize=10)
            ax.set_ylabel(r'\textbf{%s}' %ylabel, fontsize=10)
            ax.set_ylim(0, 2)
            y_ticklabels = np.arange(0, 2.5, 0.5)
            ax.set_yticks(y_ticklabels)
            ax.yaxis.set_tick_params(width=2)
            y_ticklabels = y_ticklabels.astype(str)
            for i in range(y_ticklabels.shape[0]):
                y_ticklabels[i] = r'\textbf{%s}'%y_ticklabels[i]
            ax.set_yticklabels(y_ticklabels)
            ax.set_xlim(0, self.num_days)
            xticks = np.arange(self.num_days+1)
            ax.set_xticks(xticks)
            ax.xaxis.set_tick_params(width=2)
            x_ticklabels = xticks.astype(str)
            for i in range(x_ticklabels.shape[0]):
                x_ticklabels[i] = r'\textbf{%s}'%x_ticklabels[i]
            ax.set_xticklabels(x_ticklabels)

        fig, axes = plt.subplots(nrows = np.ceil(num_concentrations * 3 /2), ncols=6)
        xlabel = "day"
        ylabel = "frequency"
        for i, ax in enumerate(axes):
            right_side = ax.spines["right"]
            top_side = ax.spines["top"]
            right_side.set_visible(False)
            top_side.set_visible(False)
            ax.set_title("Well {}".format(i), fontsize=12)
            ax.set_xlabel(r'\textbf{%s}' %xlabel, fontsize=10)
            ax.set_ylabel(r'\textbf{%s}' %ylabel, fontsize=10)
            ax.set_ylim(0, 1)
            y_ticklabels = np.arange(0, 1.2, 0.2)
            ax.set_yticks(y_ticklabels)
            ax.yaxis.set_tick_params(width=2)
            y_ticklabels = y_ticklabels.astype(str)
            for i in range(y_ticklabels.shape[0]):
                y_ticklabels[i] = r'\textbf{%s}'%y_ticklabels[i]
            ax.set_yticklabels(y_ticklabels)
            ax.set_xlim(0, self.num_days)
            xticks = np.arange(self.num_days+1)
            ax.set_xticks(xticks)
            ax.xaxis.set_tick_params(width=2)
            x_ticklabels = xticks.astype(str)
            for i in range(x_ticklabels.shape[0]):
                x_ticklabels[i] = r'\textbf{%s}'%x_ticklabels[i]
            ax.set_xticklabels(x_ticklabels)


    def generate_sequences(self, PWM_params, number_to_generate):
        '''input PWM_params of shape (num_conc, 4, num_days+1)
                 and a matched number_to_generate of shape (num_conc, )
                 seed value for np.random.seed() before running np.random.choice()
           will generate num_conc sets of sequences where each sequence has a length of num_days+1 and is built from an alphabet of (0, 1, 2, 3)
           the number of sequences in each conc set is provided by number_to_generate and is the minimum number of worms that should have been observed every day
           '''
        ## need to track it not only by concentration, but by well too. UGH

        np.random.seed(self.random_seed)
        sequences = np.full((self.num_concentrations*3, np.amax(number_to_generate), self.num_days+1, self.num_experiments), np.nan)
        to_gen_from = np.array([0, 1, 2, 3])
        for exp_index in np.arange(self.num_experiments):
            for well_index in np.arange(self.num_concentrations * 3):
                num_to_generate = number_to_generate[well_index, exp_index]
                for day_index in np.arange(self.num_days + 1): #days
                    gen_score = np.random.choice(to_gen_from, size=num_to_generate, p=PWM_params[well_index, : ,day_index, exp_index])
                    sequences[well_index, :num_to_generate, day_index, exp_index] = np.sort(gen_score)[::-1]

        return sequences

    def recreate_scores3(self, sequences):
        #need to remake self.scores3_by_conc and self.scores3_by_well
        #sequences is num_wells, num_worms_max, num_days + 1, num_experiments
        new_scores3_by_well = np.zeros((self.num_concentrations*3, 4, self.num_days, self.num_experiments)) #num_concentrations*3 because concentrations for each experiment should be in triplicate, 4 because of 0-1-2-3 scoring, num_days, and num_experiments
        new_scores3_by_well[:, 0, :, :] = np.sum(sequences[:, :, 1:, :] == 0, axis=1)
        new_scores3_by_well[:, 1, :, :] = np.sum(sequences[:, :, 1:, :] == 1, axis=1)
        new_scores3_by_well[:, 2, :, :] = np.sum(sequences[:, :, 1:, :] == 2, axis=1)
        new_scores3_by_well[:, 3, :, :] = np.sum(sequences[:, :, 1:, :] == 3, axis=1)

        new_scores3_by_conc = np.zeros((self.num_concentrations, 4, self.num_days, self.num_experiments))
        for well_index in range(new_scores3_by_well.shape[0]):
            conc_index = self.conc_index[self.well_index_to_conc[well_index]]
            new_scores3_by_conc[conc_index, :, :, :] += new_scores3_by_well[well_index, :, :, :]

        return new_scores3_by_well, new_scores3_by_conc
