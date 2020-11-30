#!/usr/bin/env python3

import numpy as np
import logging
import datetime
import pandas as pd
from scipy import stats


class WormAnalysis_probModel():
	def __init__(self, drug, strain, stage, num_days, num_experiments, num_concentrations, scores3_by_well, scores3_by_conc, total_nums, total_num_in_well, new_log_file):
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

		logging.basicConfig(filename=new_log_file, level=logging.INFO, filemode='w', format='%(name)s - %(levelname)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')

	def run(self):
		return 0

	def get_num_treated(self):
		return np.array(self.df_tablec.loc[:,r'\textbf{Number of Nematodes Treated}'])

	def learn_PWM_params(self):
		'''Here we'll use the number that are actually observed at each day to learn the params;
		when we generate based on these parameters later, we'll generate the minimum number that should have been observed each day'''
		logging.info('Because of incomplete data, learning a probabilistic model to rerun all of the analysis with generated sequences')
		logging.info('learning PWM parameters')
		params = np.zeros((self.num_concentrations, 4, self.num_days+1, self.num_experiments))
		params[:,:,0] = np.tile([0, 0, 0, 1], (self.num_concentrations, 1))#on day 0, only scores of 3

		'''use self.scores3_by_conc[conc_i, score, day_index, exp_num] to find numerator and denominator'''
		for j in range(self.num_experiments): ###NEED TO REVISIT THIS
			totals = np.sum(np.sum(self.scores3_by_conc, axis=3), axis=1).reshape((self.num_concentrations, 1, self.num_days))
			num_of_score_by_day_and_conc = np.sum(self.scores3_by_conc, axis=3)
			params[:, :, 1:] = num_of_score_by_day_and_conc/totals
		return params

	def generate_sequences(self, PWM_params, number_to_generate, seed):
		'''input PWM_params of shape (num_conc, 4, num_days+1)
				 and a matched number_to_generate of shape (num_conc, )
				 seed value for np.random.seed() before running np.random.choice()
		   will generate num_conc sets of sequences where each sequence has a length of num_days+1 and is built from an alphabet of (0, 1, 2, 3)
		   the number of sequences in each conc set is provided by number_to_generate and is the minimum number of worms that should have been observed every day
		   '''
		## need to track it not only by concentration, but by well too. UGH

		np.random.seed(seed)
		sequences = np.full((self.num_concentrations, np.amax(number_to_generate), self.num_days+1), np.nan)
		split_by_well = np.full((self.num_concentrations*3, np.amax(self.total_num_in_well), self.num_days+1, self.num_experiments), np.nan)
		to_gen_from = np.array([0, 1, 2, 3])
		for i, num in enumerate(number_to_generate):
			for j in range(PWM_params.shape[2]):
				gen_score = np.random.choice(to_gen_from, size=num, p=PWM_params[i,:,j])
				sequences[i, :num, j] = gen_score
			#randomly assign the generated paths to wells
			#use conc_to_well_index to find number that we need for each well
			num_by_well = self.total_num_in_well[self.conc_to_well_index[self.uniq_conc[i]],:]
			indexes_to_assign = np.arange(num)
			probs = np.full((indexes_to_assign.shape[0]), 1/indexes_to_assign.shape[0])
			for k in range(num_by_well.shape[0]):
				for l in range(num_by_well.shape[1]):
					assigning_indexes = np.random.choice(indexes_to_assign, size = num_by_well[k, l], replace=False, p = probs)
					split_by_well[self.conc_to_well_index[self.uniq_conc[i]][k], :num_by_well[k, l], :, l] = sequences[i, assigning_indexes, :]

					#IndexError: index 24 is out of bounds for axis 0 with size 5
					#remove assigning_indexes from indexes_to_assign by setting their probability to 0 and adjusting the rest to add to 1 again
					probs[assigning_indexes] = 0
					locs_notZero = np.where(probs != 0)[0]
					probs[locs_notZero] = 1
					probs /= locs_notZero.shape[0]


		return sequences, split_by_well

	def recreate_scores3(self, sequences_aggregated, sequences_split):
		#need to remake self.scores3_by_conc and self.scores3_by_well
		self.scores3_by_well = np.zeros((self.num_concentrations*3, 4, self.num_days, self.num_experiments)) #num_concentrations*3 because concentrations for each experiment should be in triplicate, 4 because of 0-1-2-3 scoring, num_days, and num_experiments
		self.scores3_by_conc = np.zeros((self.num_concentrations, 4, self.num_days, self.num_experiments))
		#problem! I forgot this had self.num_experiments - maybe revisit how I generated the sequence paths
