#!/usr/bin/env python3

'''
worm analysis package functions related to reporting tables
'''

def find_numTotal(consistent_num_arr, scores_by_conc, scores_by_well, nscores_insystem, num_experiments, num_concentrations, num_replicates):
    if np.all(consistent_num_arr):
        total_nums = np.sum(scores_by_conc[:,:,0,:].reshape(-1, nscores_insystem*num_experiments), axis=1).astype(np.int32)
    else:
        logging.info('Note, since there were not a consistent number of worms in a given well for every day of the experiment, we find the total number of nematodes treated by the following steps'
            + '\nfirst, finding the day with the highest number of worms in each well for each experiment'
            + '\nsecond, summing the different experiments for each well'
            + '\nfinally, summing the triplicate wells for each concentration')
        total_nums = np.sum(np.sum(np.amax(np.sum(scores_by_well, axis=1), axis=1), axis=1).reshape((num_concentrations, num_replicates)), axis=1).astype(np.int32)
    #total_num_in_well = np.amax(np.sum(scores3_by_well, axis=1), axis=1).astype(np.int32) #this is num_wells * num_experiments
    return(total_nums)
