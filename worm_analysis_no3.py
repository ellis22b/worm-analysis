#!/usr/bin/env python3

'''
worm analysis package no3 scores related functions
'''

import numpy as np

def transform_raw_to_no3(raw, nscore_insystem, num_days, num_experiments):
    transformed = np.hstack((raw[:,:nscore_insystem-1,:,:], np.sum(raw[:,nscore_insystem-1:, :, :], axis=1).reshape((-1, 1, num_days, num_experiments))))
    return (transformed)
