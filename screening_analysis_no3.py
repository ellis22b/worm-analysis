#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import logging
import datetime
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

class WormAnalysis_no3():
	def __init__(self, drug, strain, stage, uniq_conc, concUnits, concUnits_dict, conc_colors, conc_markers, conc_outline, mM, num_days, num_experiments, num_concentrations, scores3_by_well, scores3_by_conc, new_log_file):
		self.drug = drug
		self.strain = strain
		self.stage = stage
        self.uniq_conc = uniq_conc
		self.concUnits = concUnits
		self.concUnits_dict = concUnits_dict
		self.conc_colors_lo_to_hi = conc_colors
		self.conc_markers_lo_to_hi = conc_markers
		self.conc_marker_outline_lo_to_hi = conc_outline
		self.mM = mM
		self.num_days = num_days
		self.num_experiments = num_experiments
		self.num_concentrations = num_concentrations
		self.scores3_by_well = scores3_by_well
		self.scores3_by_conc = scores3_by_conc

		logging.basicConfig(filename=new_log_file, level=logging.INFO, filemode='w', format='%(name)s - %(levelname)s - %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
		

	def run(self, plotLine3, plotIT50, plotIC50, isep, expNames, rep, C_day, x0_val):
		self.transform_raw_to_no3()
		logging.info('Transformed data by combining all scores of 2&3 for each day, concentration, and experiment')
		self.motility_no3()
		logging.info('Found the motility index scores')
        self.make_tables()

		if plotLine3:
			if isep:
				for i, exp in enumerate(expNames):
					self.plotLineCurves(self.motility_index_scores_by_conc_no3[:, :, i], 'isep_motility_no3_{}_{}_{}_{}.png'.format(self.drug, self.stage, self.strain, exp), r"\textbf{%s on %s %s}" % (exp, self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$", "Motility Index Score", 0, 2, 1)
			reshaped_mid_by_well = self.motility_index_scores_by_well_no3.reshape((self.num_concentrations, 3, self.num_days+1, self.num_experiments))
			motility_index_across_exp = np.zeros((self.num_concentrations, 3*self.num_experiments, self.num_days+1), dtype=np.float64)
			for j in range(self.num_experiments):
				motility_index_across_exp[:, j*3:(j*3)+3, :] = reshaped_mid_by_well[:, :, :, j]
			motility_avg_across_exp = np.mean(motility_index_across_exp, axis=1)
			self.plotLineCurves(motility_avg_across_exp, 'average_mortality_no3_{}_{}_{}.png'.format(self.drug, self.stage, self.strain), r"\textbf{%s on %s %s}" %(self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$", "Motility Index Score", 0, 2, 1)

		if plotIT50:
			logging.info('Beginning Survival Analysis for {} as the representative experiment'.format(expNames[rep-1]))
            inhibited = np.full((self.num_concentrations, self.num_days+1), 100, dtype=np.float64)
            toAnalyze = self.scores2_by_conc[:, :, :, rep-1]
            num_total = np.sum(toAnalyze[:, :, 0], axis=1).reshape((self.num_concentrations, 1))
            num_at_risk = to_Analyze[:, 2,:] #num at risk is only worms of score 2
            num_at_risk_corrected = np.minimum.accumulate(num_at_risk, axis=1)
            inhibited[:, 1:] = num_at_risk_corrected/num_total*100
            figname_base = '{}_IT50' + '{}_{}_{}_{}.png'.format(self.drug, self.stage, self.strain, expNames[rep-1])
            self.plotSurvivalTime(inhibited, 0, figname_base, plot_mortality=False)
            T50 = np.sum(toPopulate <= 50.0, axis=1)
            T50 = (self.num_days + 1 - T50).astype(str)
            bool_U = T50 == str(self.num_days + 1)
            T50[bool_U] = 'U'
            self.df_tablec[r'\textbf{IT50}'] = T50
            logging.info("Added the IT50 values to the table which will be printed later. Column name is 'IT50'")

        if plotIC50:
            print('not done yet')

    def make_tables():
        index = []
        units = np.tile(self.concUnits_dict[self.concUnits], self.uniq_conc.shape[0])
        for conc, unit in zip(self.uniq_conc, units):
            index.append(r'\textbf{%s %s}' % (conc,unit))
        self.df_tablec = pd.DataFrame(index=index)
        self.df_tablec.index.name='Concentration'

        days_arr = np.arange(self.num_days+1).astype(str)
        for k in range(days_arr.shape[0]):
            days_arr[k] = r'\textbf{%s}' % days_arr[k]
        self.df_tabled = pd.DataFrame(index=days_arr)
        self.df_tabled.index.name='Day'
        self.df_tabled[r'\textbf{IC50}'] = np.tile('NC', days_arr.shape[0])

	def transform_raw_to_no3(self):
		self.scores2_by_well = np.hstack((self.scores3_by_well[:, :2, :, :], np.sum(self.scores3_by_well[:, 2:, :, :], axis=1).reshape((-1, 1, self.num_days, self.num_experiments))))
		self.scores2_by_conc = np.hstack((self.scores3_by_conc[:, :2, :, :], np.sum(self.scores3_by_conc[:, 2:, :, :], axis=1).reshape((self.num_concentrations, 1, self.num_days, self.num_experiments))))

	def motility_no3(self):
		self.motility_index_scores_by_conc_no3 = np.full((self.num_concentrations, self.num_days+1, self.num_experiments), 2.0, dtype=np.float64) #day 0 will be set to 2 automatically, will set days1 onwards at the end
		adjusted_sum = np.sum(self.scores2_by_conc * np.array([0, 1, 2]).reshape((1,-1,1,1)), axis=1) #weight number of worms by the score and sum: 0*num0_in_well + 1*num1_in_well + 2*num2|3_in_well
		divided = adjusted_sum/np.sum(self.scores2_by_conc, axis=1) #divide the sum above (adjusted_sum) by the number of worms total in the well (the denominator)
		self.motility_index_scores_by_conc_no3[:, 1:, :] = divided #setting days1 onwards

		self.motility_index_scores_by_well_no3 = np.full((self.num_concentrations*3, self.num_days+1, self.num_experiments), 2.0, dtype=np.float64)
		adjusted_sum_by_well = np.sum(self.scores2_by_well * np.array([0, 1, 2]).reshape((1,-1,1,1)), axis=1)
		divided_by_well = adjusted_sum_by_well/np.sum(self.scores2_by_well, axis=1)
		self.motility_index_scores_by_well_no3[:, 1:, :] = divided_by_well

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