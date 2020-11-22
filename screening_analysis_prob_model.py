#!/usr/bin/env python3

import numpy as np
import logging
import datetime
import pandas as pd
from scipy import stats
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