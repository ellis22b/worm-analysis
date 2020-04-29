## Worm Analysis
This repository contains the code for analysis of health-rating *C. elegans* drug screening experiments. All code is based off the analysis used in [this paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0179376)<br />

### Required inputs
Use the `--toAssess` argument to pass the experiment input files. (`--toAssess kw9_export.txt kw11_export.txt kw12_export.txt` will pass 3 input files). <br /><br />

The input files should be tab delimited and formatted like this example beginning 37 lines of `kw9_export.txt` where each line specifies a well at a given concentration on a given day and the number of worms that were scored as 0, 1, 2, or 3 for that well/day/concentration as well as the total number of worms in the well, and a boolean argument for the final column in only the first day's set of lines.
```
Well	Concentration	Day	num0_in_well	num1_in_well	num2_in_well	num3_in_well	numTotal_in_well	numTotal_equal
1	0	1	0	0	0	8	8	TRUE
2	0	1	0	0	0	7	7	TRUE
3	0	1	0	0	0	6	6	TRUE
4	0.1	1	0	0	0	7	7	TRUE
5	0.1	1	0	0	1	6	7	TRUE
6	0.1	1	0	0	0	4	4	TRUE
7	1	1	0	0	0	8	8	TRUE
8	1	1	0	0	0	5	5	TRUE
9	1	1	0	0	0	6	6	TRUE
10	10	1	0	0	2	4	6	TRUE
11	10	1	0	0	0	10	10	TRUE
12	10	1	0	0	1	6	7	TRUE
13	100	1	0	0	0	3	3	TRUE
14	100	1	0	0	0	8	8	TRUE
15	100	1	0	0	0	12	12	TRUE
16	1000	1	0	0	4	1	5	TRUE
17	1000	1	0	0	0	12	12	TRUE
18	1000	1	0	0	0	6	6	TRUE
1	0	2	0	0	0	8	8
2	0	2	0	0	0	7	7
3	0	2	0	0	0	6	6
4	0.1	2	0	0	0	7	7
5	0.1	2	0	0	0	7	7
6	0.1	2	0	0	0	4	4
7	1	2	0	0	0	8	8
8	1	2	0	0	0	5	5
9	1	2	0	0	0	6	6
10	10	2	0	0	3	3	6
11	10	2	0	0	0	10	10
12	10	2	1	0	1	5	7
13	100	2	0	0	3	0	3
14	100	2	0	0	2	6	8
15	100	2	0	0	0	12	12
16	1000	2	0	0	4	1	5
17	1000	2	0	0	0	12	12
18	1000	2	0	0	5	1	6
```

If performing IT50/LT50 analysis, `--representative` is a **required** argument.

### Note on default analysis
Note that the default is to run and plot all analysis including the following. In order to not assess any of the following, set the specified argument to false ( `--plotLine3 False` will cause the program to skip plotting motility index score)
- motility index score by concentration over the length of the experiment (`--plotLine3`)
- mortality by concentration over the length of the experiment (`--plotLine1`)
- IC50 (`--plotIC50`)
- LC50 (`--plotLC50`)
- IT50 (`--plotIT50`)
- LT50 (`--plotLT50`)

### IC50/LC50 specific notes and arguments
If running IC50 and/or LC50, the default is to consider day 4 data. The user can change this using the `--C_day` argument. `--C_day 5` would analyze day 5. <br /><br />
For IC50 and LC50 experiments, the process is based off of the [prism equation](https://www.graphpad.com/guides/prism/8/curve-fitting/reg_dr_inhibit_variable.htm) for log(inhibitor) vs response -- variable slope. Concentration is log10 transformed replacing X=[0] with a default of X=[1e-6] or a user defined value using the argument `--x0_value`.

### IT50/LT50 specific notes and required argument
If running, IT50 and/or LT50, the user must define the representative experiment to run, using the `--representative` argument. Specifically, provide the number corresponding to the index (1 based index) of the representative experiment in the `--toAssess` argument. If `--toAssess kw9_export.txt kw11_export.txt kw12_export.txt` and the user wants kw12 to be the representative experiment, then use `--representative 3`. For kw11 as representative use `--representative 2`. If a representative experiment is not defined by the user, the script will exit with a message after performing all other analysis. <br /><br />

### General notes on the program
The number of days and concentrations of a drug in an experiment are inferred from the input data. Any inconsistencies in a well number's concentration across experiments raises an error. <br /><br />

A logfile stores the user defined settings to run the analysis, warnings, and info pertaining to derived values and parameters from curve-fitting. The logfile name can be specified using the `--logfile` argument. Otherwise the default is set to reflect the date and time that the analysis was run.<br /><br />
