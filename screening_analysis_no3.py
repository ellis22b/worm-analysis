#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import sem
from scipy.optimize import curve_fit
from scipy.interpolate import make_interp_spline, BSpline
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
	def __init__(self, consistent_num, drug, strain, stage, uniq_conc, concUnits, concUnits_dict, conc_index, well_index_to_conc, conc_colors, conc_markers, conc_outline, C_day, x0_val, mM, num_days, num_experiments, num_concentrations, scores3_by_well, scores3_by_conc):
		self.consistent_num = consistent_num
		self.drug = drug
		self.strain = strain
		self.stage = stage
		self.uniq_conc = uniq_conc
		self.concUnits = concUnits
		self.concUnits_dict = concUnits_dict
		self.conc_index = conc_index
		self.well_index_to_conc = well_index_to_conc
		self.conc_colors_lo_to_hi = conc_colors
		self.conc_markers_lo_to_hi = conc_markers
		self.conc_marker_outline_lo_to_hi = conc_outline
		self.C_day = C_day
		self.x0_val = x0_val
		self.mM = mM
		self.num_days = num_days
		self.num_experiments = num_experiments
		self.num_concentrations = num_concentrations
		self.scores3_by_well = scores3_by_well
		self.scores3_by_conc = scores3_by_conc

	def run(self, plotLine3, plotIT50, plotIC50, isep, expNames, rep, C_day, x0_val, hill2, spline_k1, spline_k2, default=True, not_default_2={}):
		self.transform_raw_to_no3()
		logging.info('Transformed data by combining all scores of 2&3 for each day, concentration, and experiment')
		self.motility_no3()
		logging.info('Found the motility index scores')
		self.make_tables()

		if plotLine3:
			if isep:
				for i, exp in enumerate(expNames):
					self.plotLineCurves(self.motility_index_scores_by_conc_no3[:, :, i], 'isep_motility_{}_{}_{}_{}_no3.png'.format(self.drug, self.stage, self.strain, exp), r"\textbf{%s on %s %s}" % (exp, self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$", "Motility Index Score", 0, 2, 1)
			reshaped_mid_by_well = self.motility_index_scores_by_well_no3.reshape((self.num_concentrations, 3, self.num_days+1, self.num_experiments))
			motility_index_across_exp = np.zeros((self.num_concentrations, 3*self.num_experiments, self.num_days+1), dtype=np.float64)
			for j in range(self.num_experiments):
				motility_index_across_exp[:, j*3:(j*3)+3, :] = reshaped_mid_by_well[:, :, :, j]
			motility_avg_across_exp = np.mean(motility_index_across_exp, axis=1)
			self.plotLineCurves(motility_avg_across_exp, 'average_motility_{}_{}_{}_no3.png'.format(self.drug, self.stage, self.strain), r"\textbf{%s on %s %s}" %(self.drug, self.stage, self.strain) + r" $\textbf{$\textit{C. elegans}$}$", "Motility Index Score", 0, 2, 1)

		if plotIT50:
			logging.info('Beginning Survival Analysis for {} as the representative experiment'.format(expNames[rep-1]))
			inhibited = np.full((self.num_concentrations, self.num_days+1), 100, dtype=np.float64)
			toAnalyze = self.scores2_by_conc[:, :, :, rep-1]
			if self.consistent_num:
				num_total = np.sum(toAnalyze[:, :, 0], axis=1).reshape((self.num_concentrations, 1))
			else:
				num_total = np.sum(toAnalyze, axis=1)
			num_at_risk = toAnalyze[:, 2,:] #num at risk is only worms of score 2
			num_at_risk_corrected = np.minimum.accumulate(num_at_risk, axis=1)
			inhibited[:, 1:] = num_at_risk_corrected/num_total*100
			figname_base = '{}_IT50' + '{}_{}_{}_{}_no3.png'.format(self.drug, self.stage, self.strain, expNames[rep-1])
			self.plotSurvivalTime(inhibited, 0, figname_base, plot_mortality=False)
			T50 = np.sum(inhibited <= 50.0, axis=1)
			T50 = (self.num_days + 1 - T50).astype(str)
			bool_U = T50 == str(self.num_days + 1)
			T50[bool_U] = 'U'
			self.df_tablec[r'\textbf{IT50}'] = T50
			logging.info("Added the IT50 values to the table which will be printed later. Column name is 'IT50'")

		if plotIC50:
			totals = np.sum(self.scores2_by_well[:,:,self.C_day-1,:], axis=1)
			uninhibited2 = self.scores2_by_well[:,2,self.C_day-1,:]/totals*100

			toAvg2 = np.zeros((self.num_concentrations, 3*self.num_experiments))
			for i in range(uninhibited2.shape[0]):
				toAvg2[self.conc_index[self.well_index_to_conc[i]], ((i%3)*self.num_experiments): (((i%3)+1) * self.num_experiments)] = uninhibited2[i].copy()

			sem2 = sem(toAvg2, axis=1)
			avg2 = np.average(toAvg2, axis=1)

			conc_X = np.tile(np.array([self.well_index_to_conc[x] for x in np.arange(self.num_concentrations*3)]).reshape(-1,1), (1, self.num_experiments))

			if self.num_experiments >= 3 and default:
				param_dict = self.set_guess_params(avg2, hill2)
			elif self.num_experiments >= 3 and not default:
				param_dict = self.set_guess_params(avg2, hill2, default=False, not_default_2 = not_default_2)
			elif self.num_experiments < 3 and default:
				param_dict = self.set_guess_params(avg2, hill2, num_exper='l3')
			elif self.num_experiments < 3 and not default:
				param_dict = self.set_guess_params(avg2, hill2, default=False, not_default_2 = not_default_2, num_exper='l3')

			P0_20, P0_20_top, P0_20_bottom, P0_20_ic50, P0_20_hill = param_dict['P0_20'], param_dict['P0_20_top'], param_dict['P0_20_bottom'], param_dict['P0_20_ic50'], param_dict['P0_20_hill']

			lowest_nz_conc = np.sort(self.uniq_conc)[1]
			highest_conc = np.amax(self.uniq_conc)

			if self.num_experiments >= 3:
				logging.info('Running Levenberg-Marquardt Algorithm Scipy Curve Fitting for 2-1-0 scoring using the default max number of function evaluations. Initial values are the following. \nINITIAL Top:\t{}\nINITIAL Bottom:\t{}\nINITIAL IC50:\t{}\nINITIAL HillSlope:\t{}'.format(P0_20_top, P0_20_bottom, P0_20_ic50, P0_20_hill))
				popt2, popc2 = curve_fit(self.inhibitorResponse_equation, conc_X.flatten(), uninhibited2.flatten(), p0=P0_20, method='lm')
				top_20, bottom_20, ic50_20, hill_20 = popt2[0], popt2[1], popt2[2], popt2[3]
				logging.info('Returned lm fit for 2-1-0 scoring.\nFIT Top:\t{}\nFIT Bottom:\t{}\nFIT IC50:\t{}\nFIT HillSlope:\t{}'.format(top_20, bottom_20, ic50_20, hill_20))
				no_fit_bool = self.evaluate_no_fit(top_20, bottom_20, ic50_20, highest_conc)
				self.set_to_log_value(bottom_20, ic50_20, highest_conc, lowest_nz_conc)
				if no_fit_bool:
					self.plotIC_noFit(r'$\mathrm{IC_{50}}$' + ' {} on {} {} Day {}'.format(self.drug, self.stage, self.strain, self.C_day), 'IC50_{}_{}_{}_{}_no3.png'.format(self.drug, self.stage, self.strain, self.C_day), self.uniq_conc, avg2, sem2, spline_k1, spline_k2)
				else:
					self.plotIC(r'$\mathrm{IC_{50}}$' + ' {} on {} {} Day {}'.format(self.drug, self.stage, self.strain, self.C_day), 'IC50_{}_{}_{}_{}_no3.png'.format(self.drug, self.stage, self.strain, self.C_day), self.uniq_conc, avg2, sem2, ic50_20, curve_fit_hillslope=hill_20, curve_fit_top=top_20, curve_fit_bottom=bottom_20)
				logging.info('Completed Non-linear Regression for Inhibition Response Analysis')
			else:
				logging.info('Running Levenberg-Marquardt Algorithm Scipy Curve Fitting for 2-1-0 scoring useing the default max number of function evaluations and a constant hill slope of {}. Initial values are the following.\nINITIAL Top:\t{}\nINITIAL Bottom:\t{}\nINITIAL IC50:\t{}'.format(P0_20_hill, P0_20_top, P0_20_bottom, P0_20_ic50))
				popt2, popc2 = curve_fit(self.inhibitorResponse_equation, conc_X.flatten(), uninhibited2.flatten(), p0=P0_20, method='lm')
				top_20, bottom_20, ic50_20 = popt2[0], popt2[1], popt2[2]
				logging.info('Returned lm fit for 2-1-0 scoring.\nFIT Top:\t{}\nFIT Bottom:\t{}\nFIT IC50:\t{}'.format(top_20, bottom_20, ic50_20))
				no_fit_bool = self.evaluate_no_fit(top_20, bottom_20, ic50_20, highest_conc)
				self.set_to_log_value(bottom_20, ic50_20, highest_conc, lowest_nz_conc)
				if no_fit_bool:
					self.plotIC_noFit(r'$\mathrm{IC_{50}}$' + ' {} on {} {} Day {}'.format(self.drug, self.stage, self.strain, self.C_day), 'IC50_{}_{}_{}_{}_no3.png'.format(self.drug, self.stage, self.strain, self.C_day), self.uniq_conc, avg2, sem2, spline_k1, spline_k2)
				else:
					self.plotIC(r'$\mathrm{IC_{50}}$' + ' {} on {} {} Day {}'.format(self.drug, self.stage, self.strain, self.C_day), 'IC50_{}_{}_{}_{}_no3.png'.format(self.drug, self.stage, self.strain, self.C_day), self.uniq_conc, avg2, sem2, ic50_20, curve_fit_hillslope=P0_20_hill, curve_fit_top=top_20, curve_fit_bottom=bottom_20)
				logging.info('Completed Non-linear Regression for Inhibition Response Analysis but with constrained hill slope and therefore the fit is likely less than ideal')
		self.reportTable(expNames[rep-1], plotIT50, plotIC50)

	def set_to_log_value(self, bottom_20, ic50_20, highest_conc, lowest_nz_conc):
		logging_value = r'\textbf{IC50}'
		logging_value2 = 'IC50'
		logging_day = r'\textbf{%s}' % str(self.C_day)
		to_log = '{0:.3f}'.format(ic50_20)
		if bottom_20 > 50 or ic50_20 > highest_conc:
			to_log = r'\textgreater' + '{}'.format(highest_conc)
		elif bottom_20 < 50 and ic50_20 < lowest_nz_conc:
			to_log = r'\textless' + '{}'.format(lowest_nz_conc)
		self.df_tabled.loc[logging_day, logging_value] = to_log
		logging.info("Added the %s value to the table which will be printed later. Column name is '%s' and the day/row is '%s" %(logging_value2, logging_value2, self.C_day))

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
		exponent = (np.log10(ic50)-self.transform_X(X))*hillslope
		equation = bottom + (top-bottom)/(1+(10**exponent))
		return equation

	def set_guess_params(self, avg2, hill2, default=True, num_exper='p3', not_default_2={}):
		if default:
			P0_20_top = np.amax(avg2)
			P0_20_bottom = np.amin(avg2)
			P0_20_ic50 = 1
			if num_exper == 'p3':
				P0_20_hill = -1
				P0_20 = [P0_20_top, P0_20_bottom, P0_20_ic50, P0_20_hill]
				return {'P0_20': P0_20,
						'P0_20_top': P0_20_top,
						'P0_20_bottom': P0_20_bottom,
						'P0_20_ic50': P0_20_ic50,
						'P0_20_hill': P0_20_hill}
			elif num_exper == 'l3':
				P0_20 = [P0_20_top, P0_20_bottom, P0_20_ic50]
				return {'P0_20': P0_20,
						'P0_20_top': P0_20_top,
						'P0_20_bottom': P0_20_bottom,
						'P0_20_ic50': P0_20_ic50,
						'P0_20_hill': hill2}
		else:
			if num_exper == 'p3':
				P0_20 = [not_default_2['top'], not_default_2['bottom'], not_default_2['ic50'], not_default_2['hill']]
				return {'P0_20': P0_20,
						'P0_20_top': not_default_2['top'],
						'P0_20_bottom': not_default_2['bottom'],
						'P0_20_ic50': not_default_2['ic50'],
						'P0_20_hill': not_default_2['hill']}
			elif num_exper == 'l3':
				P0_20 = [not_default_2['top'], not_default_2['bottom'], not_default_2['ic50']]
				return {'P0_20': P0_20,
						'P0_20_top': not_default_2['top'],
						'P0_20_bottom': not_default_2['bottom'],
						'P0_20_ic50': not_default_2['ic50'],
						'P0_20_hill': hill2}

	def make_tables(self):
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

	def dose_response_sigmoid(self, X, X_mid, hill, bottom, top):
		Y = (bottom+(top-bottom)/(1+np.exp(hill*(np.log(X)-np.log(X_mid)))))
		return Y

	def plotIC(self, title, figname, concs_toAnalyze, averages, sems, curve_fit_ic50, curve_fit_hillslope=-1, curve_fit_top=100, curve_fit_bottom=0, ylabel='\% Uninhibited', ysep=20, ymin=0, ymax=100):
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

		curve1e = self.dose_response_sigmoid(linspace_x1_antilog, curve_fit_ic50, curve_fit_hillslope, curve_fit_top, curve_fit_bottom)
		curve2e = self.dose_response_sigmoid(linspace_x2_antilog, curve_fit_ic50, curve_fit_hillslope, curve_fit_top, curve_fit_bottom)

		fig = plt.figure(constrained_layout=False)
		widths=[1,8]
		gs = fig.add_gridspec(1,2, width_ratios=widths, wspace=0.05)
		ax1 = fig.add_subplot(gs[0,0])
		ax2 = fig.add_subplot(gs[0,1])

		ax1.axhline(50, linestyle=':', color='black', clip_on=False)
		ax2.axhline(50, linestyle=':', color='black', clip_on=False)

		for side in ['right', 'top']:
			spine_side1 = ax1.spines[side]
			spine_side2 = ax2.spines[side]
			spine_side1.set_visible(False)
			spine_side2.set_visible(False)
		left2 = ax2.spines['left']
		left2.set_visible(False)

		fig.suptitle(r'\textbf{%s}' %title)
		ax1.set_ylabel(r'\textbf{%s}' %ylabel)
		ax1.set_ylim(ymin, ymax, ysep)
		ax2.set_ylim(ymin, ymax, ysep)
		ax2.set_yticks([])
		ax2.set_yticklabels([])
		xlabel='Concentration (%s)' % self.concUnits_dict[self.concUnits]
		fig.text(0.5, 0.02, r'\textbf{%s}' %xlabel, ha='center', va='center')

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

		ax1.plot((1,1), (1,1), color='black', clip_on=False)
		ax2.plot((0,0), (1,1), color='black', clip_on=False)

		fig.savefig(figname)
		plt.close(fig)
		logging.info('Plotted the figure {}'.format(figname))

	def plotIC_noFit(self, title, figname, concs_toAnalyze, averages, sems, spline_k1, spline_k2, ylabel='\% Uninhibited', ysep=20, ymin=0, ymax=100):
		conc_ticks = np.copy(self.uniq_conc)
		bool_0 = conc_ticks == 0
		conc_ticks[bool_0] = self.x0_val
		log_conc_ticks = np.log10(conc_ticks)

		concs = np.copy(concs_toAnalyze)
		bool_0 = concs == 0
		concs[bool_0] = self.x0_val
		log_concs = np.log10(concs)

		fig = plt.figure(constrained_layout=False)
		widths=[1,8]
		gs = fig.add_gridspec(1,2,width_ratios=widths, wspace=0.05)
		ax1 = fig.add_subplot(gs[0,0])
		ax2 = fig.add_subplot(gs[0,1])

		ax1.axhline(50, linestyle=':', color='black', clip_on=False)
		ax2.axhline(50, linestyle=':', color='black', clip_on=False)

		for side in ['right', 'top']:
			spine_side1 = ax1.spines[side]
			spine_side2 = ax2.spines[side]
			spine_side1.set_visible(False)
			spine_side2.set_visible(False)
		left2 = ax2.spines['left']
		left2.set_visible(False)

		fig.suptitle(r'\textbf{%s}' %title)
		ax1.set_ylabel(r'\textbf{%s}' %ylabel)
		ax1.set_ylim(ymin, ymax, ysep)
		ax2.set_ylim(ymin, ymax, ysep)
		ax2.set_yticks([])
		ax2.set_yticklabels([])
		xlabel = 'Concentration (%s)' % self.concUnits_dict[self.concUnits]
		fig.text(0.5, 0.02, r'\textbf{%s}' %xlabel, ha='center', va='center')

		e = ax1.errorbar(log_concs[0], averages[0], yerr=sems[0], marker='o', fmt='none', linestyle='None', color='black', capsize=5, clip_on=False)
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
		s = ax2.scatter(log_conc_ticks[1:], averages[1:], c='black')
		ax2.plot(x_s_2, power_smooth_2, c='black', clip_on=False)

		s.set_clip_on(False)

		ax2.set_xlim(log_conc_ticks[1], log_conc_ticks[-1])
		ax2.set_xticks(log_conc_ticks[1:])
		ticklabels = np.hstack((conc_ticks[1], conc_ticks[2:].astype(np.int32).astype(str)))
		ax2.set_xticklabels(ticklabels)

		ax1.plot((1,1), (1,1), color='black', clip_on=False)#bottom-left line
		ax2.plot((0,0), (1,1), color='black', clip_on=False)#bottom-right line

		fig.savefig(figname)
		plt.close(fig)
		logging.info('Plotted the figure {}'.format(figname))


	def reportTable(self, rep_exp, plotIT50, plotIC50, figname_base='_no3'):
		if plotIT50 and plotIC50:
			fig, (ax1, ax2) = plt.subplots(2,1)
			#hide axes
			fig.patch.set_visible(False)
			ax1.axis('off')
			ax2.axis('off')
			ax1.axis('tight')
			ax2.axis('tight')
			#write table
			ax1.table(cellText=self.df_tablec.values, colLabels=self.df_tablec.columns, rowLabels=self.df_tablec.index, cellLoc='center', loc='center')
			ax2.table(cellText=self.df_tabled.values, colLabels=self.df_tabled.columns, rowLabels=self.df_tabled.index, cellLoc='center', loc='center')
			filename='table_{}_{}_{}_t50repexp_{}_c50day_{}{}.pdf'.format(self.drug, self.stage, self.strain, rep_exp, self.C_day, figname_base)

		elif plotIT50 and not plotIC50:
			fig, ax = plt.subplots()
			fig.patch.set_visible(False)
			ax.axis('off')
			ax.axis('tight')
			ax.table(cellText=self.df_tablec.values, colLabels=self.df_tablec.columns, rowLabels=self.df_tablec.index, cellLoc='center', loc='center')
			filename='table_{}_{}_{}_t50repexp_{}{}.pdf'.format(self.drug, self.stage, self.strain, rep_exp, figname_base)

		elif plotIC50 and not plotIT50:
			fig, ax = plt.subplots()
			fig.patch.set_visible(False)
			ax.axis('off')
			ax.axis('tight')
			ax.table(cellText=self.df_tabled.values, colLabels=self.df_tabled.columns, rowLabels=self.df_tabled.index, cellLoc='center', loc='center')
			filename='table_{}_{}_{}_c50day_{}{}.pdf'.format(self.drug, self.stage, self.strain, self.C_day, figname_base)

		fig.tight_layout()
		fig.savefig(filename, format='pdf')
		plt.close(fig)
		logging.info('Wrote table(s) of computed values to {}'.format(filename))
