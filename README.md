## Worm Analysis
This repository contains the code for analysis of health-rating *C. elegans* drug screening experiments. All code is based off the analysis used in [this paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0179376)<br />

### Required inputs
Use the `--toAssess` argument to pass the experiment/data input files. (`--toAssess kw9_export.txt kw11_export.txt kw12_export.txt` will pass 3 input files). <br /><br />

If performing IT50/LT50 analysis, `--representative` is a **required** argument. (ex: `--representative 3`). See [this part of the documentation](https://github.com/ellis22b/worm-analysis#it50lt50-specific-notes-and-required-argument) for specifics and explanation about this argument.

The following set of inputs will be used to annotate plots from the data passed to `--toAssess`
- Use `--strain` to specify the strain of *C. elegans* was treated in the assays
- Use `--lifestage` to specify which stage of *C. elegans* was treated in the assays
- Use `--drug` to specify which drug was used for treatment in the assays (concentration is assumed to be ug/mL)
- Use `--expNames` to specify the names of the experiments to be assessed (ex: `--expNames kw9 kw11 kw12`). Order of the names should match the order of the files.
- Use `--concUnits` to specify the units of the concentrations. The actual values will be inferred from the input files. `--concUnits 0` would use ug/mL in labeling plots and tables. `--concUnits 7` would use mg/kg in labeling plots and tables
  - 0: ug/mL
  - 1: uM
  - 2: M
  - 3: mM
  - 4: nM
  - 5: ng/mL
  - 6: mg/mL
  - 7: mg/kg
  - 8: mg/g

  
Another important note about the data input files is that **all inputs (at this time) in a single run should be experiments with the same strain, life stage, drug, number of days, and concentrations**. For example, in the `--toAssess` list above, kw9, kw11, and kw12 were N2 L4 Albendazole 7-day experiments with concentrations 0, 0.1, 1, 10, 100, 1000 ug/mL. (Later additions to the program can explore ways to plot drug screening more efficiently)<br />< br />

Concerning **format for the input files**: These input files should be **tab delimited** and formatted like this example beginning 37 lines of `kw9_export.txt` where each line specifies a well at a given concentration on a given day and the number of worms that were scored as 0, 1, 2, or 3 for that well/day/concentration as well as the total number of worms in the well, and a boolean argument for the final column in only the first day's set of lines.
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

### Note on default analysis
Note that the default is to run and plot all analysis including the following. In order to not assess any of the following, set the specified argument to false ( `--plotLine3 False` will cause the program to skip plotting motility index score)
- motility index score by concentration over the length of the experiment (`--plotLine3`)
- mortality by concentration over the length of the experiment (`--plotLine1`)
- IC50 (`--plotIC50`)
- LC50 (`--plotLC50`)
- IT50 (`--plotIT50`)
- LT50 (`--plotLT50`)
- report the number of worms (`--reportNum`)

### Line Plots of daily motility and mortality response by concentration specific notes and arguments
In the [2017 paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0179376), these graphs combined the scores for each individual well in each experiment for a given concentration to find the average score. These are the plots that will be produced by default by this program if `--plotLine3` and/or `--plotLine1` are true.<br /><br />
In addition to this default behavior, the argument `--include_single_exp_plots` can be used to produce plots for each individual experiment. `--include_single_exp_plots` is false by default. When true (in combination with `--plotLine3` and/or `--plotLine1`), the program will also plot each experiment by itself in addition to the default plots.

### IC50/LC50 specific notes and arguments
If running IC50 and/or LC50, the default is to consider day 4 data. The user can change this using the `--C_day` argument. `--C_day 5` would analyze day 5. <br /><br />
For IC50 and LC50 experiments, the process is based off of the [prism equation](https://www.graphpad.com/guides/prism/8/curve-fitting/reg_dr_inhibit_variable.htm) for log(inhibitor) vs response -- variable slope. Concentration is log10 transformed replacing X=[0] with a default of X=[1e-6] or a user defined value using the argument `--x0_value`.

### IT50/LT50 specific notes and required argument
If running, IT50 and/or LT50, the user must define the representative experiment to run, using the `--representative` argument. Specifically, provide the number corresponding to the index (1 based index) of the representative experiment in the `--toAssess` argument. If `--toAssess kw9_export.txt kw11_export.txt kw12_export.txt` and the user wants kw12 to be the representative experiment, then use `--representative 3`. For kw11 as representative use `--representative 2`. If a representative experiment is not defined by the user, the script will exit with a message after performing all other analysis. <br /><br />

### General notes on the program
The argument `--reportNum` is true by default and will report the total number of worms assessed by concentration (corresponding to [Table 1 of the 2017 paper](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0179376)). These total numbers will be reported in a text file named `table1_equiv_numTreated_[drug]_[stage]_[strain].txt`, where drug, stage, and strain are dependent on the inputs. Specifically if `table1_equiv_numTreated_ALB_L4_N2.txt` would be the file from input arguments `--drug ALB --lifestage L4 --strain N2`. The logfile will also record the name of the file written too.

The number of days and concentrations of a drug in an experiment are inferred from the input data. Any inconsistencies in a well number's concentration across experiments raises an error. <br /><br />

A logfile stores the user defined settings to run the analysis, warnings, and info pertaining to derived values and parameters from curve-fitting. The logfile name can be specified using the `--logfile` argument. Otherwise the default is set to reflect the date and time that the analysis was run.<br /><br />
