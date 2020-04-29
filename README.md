## Worm Analysis
This repository contains the code for analysis of health-rating *C. elegans* drug screening experiments. All code is based off the analysis used in [this paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0179376)<br />

Use the `--toAssess` argument to pass the experiment input files. (`--toAssess kw9_export.txt kw11_export.txt kw12_export.txt` will pass 3 input files). <br /><br />

Note that the default is to run and plot all analysis including the following. In order to not assess any of the following, set the specified argument to false ( `--plotLine3 False` will cause the program to skip plotting motility index score)
- motility index score by concentration over the length of the experiment (`--plotLine3`)
- mortality by concentration over the length of the experiment (`--plotLine1`)
- IC50 (`plotIC50`)
- LC50 (`plotLC50`)
- IT50 (`plotIT50`)
- LT50 (`plotLT50`)

If running IC50 and/or LC50, the default is to consider day 4 data. The user can change this using the `--C_day` argument. `--C_day 5` would analyze day 5. <br /><br />

If running, IT50 and/or LT50, the user must define the representative experiment to run, using the `--representative` argument. Specifically, provide the number corresponding to the index (1 based index) of the representative experiment in the `--toAssess` argument. If `--toAssess kw9_export.txt kw11_export.txt kw12_export.txt` and the user wants kw12 to be the representative experiment, then use `--representative 3`. For kw11 as representative use `--representative 2`. If a representative experiment is not defined by the user, the script will exit with a message after performing all other analysis. <br /><br />

The number of days and concentrations of a drug in an experiment are inferred from the input data. Any inconsistencies in a well number's concentration across experiments raises an error. <br /><br />

A logfile stores the user defined settings to run the analysis, warnings, and info pertaining to derived values and parameters from curve-fitting. The logfile name can be specified using the `--logfile` argument. Otherwise the default is set to reflect the date and time that the analysis was run.

For IC50 and LC50 experiments, the process is based off of the [prism equation](https://www.graphpad.com/guides/prism/8/curve-fitting/reg_dr_inhibit_variable.htm) for log(inhibitor) vs response -- variable slope. Concentration is log10 transformed replacing X=[0] with a default of X=[1e-6] or a user defined value using the argument `--x0_value`.
