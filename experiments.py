"""Run experiments over a range of parameters of interest"""

# TODO: run some parameter sweeps on: 
    # * neighborhood CC (in relationship to tree cc)
    # * growth rate
    # * carrying capacity
    # * ignition probability
    # collect data on mean and sd of biomass and patch size (from 200-100 years or something)
import numpy as np
from simulation import run_simulation
import pandas as pd
# TEST final experiments 
# neighborhood cc

# DEFAULT PARAMS
m = 0.5 # placeholder, not used. 
L = 30
time_steps = 200
d = 4
init_grass=0.2 
init_tree=0.1
p_disp=.05
p_prop=0.1
min_seed=2
r_grow=0.8
tree_carrying_capacity = 50
neighborhood_carrying_capacity = 100
max_ignite=0.01
rand_seed = 1551

num_reps = 2
burnout_time = 100

rng = np.random.default_rng(rand_seed)

# Neighborhood CC 
# Run over a range of different ratios of tree_carrying_capacity/neighborhood_carrying_capacity
# from 1/10 to 1 (1/10, 1/9, 1/8, etc)

# set up data collection
# will have a couple of summary stats for each iteration

ratios = range(20)
num_trials = len(ratios) * num_reps
ratio = np.zeros(num_trials)
mean_biomass = np.zeros(num_trials)
median_biomass = np.zeros(num_trials)
max_biomass = np.zeros(num_trials)
min_biomass = np.zeros(num_trials)
std_biomass = np.zeros(num_trials)
mean_patch_count = np.zeros(num_trials)
mean_mean_patch_size = np.zeros(num_trials)
std_mean_patch_size = np.zeros(num_trials)
ln_num_avg_patch_int = np.zeros(num_trials)
ln_num_avg_patch_slope = np.zeros(num_trials)

# keep track of trial number
trial_counter = 0

for r in ratios: 
    r = r + 1
    n_cc = tree_carrying_capacity*r
    for rep in range(num_reps):
        print(trial_counter)
        results_dict = run_simulation(m=m, 
                L=L, 
                t_steps=time_steps, 
                d=d, 
                init_grass=init_grass, 
                init_tree=init_tree, 
                p_disp=p_disp, 
                p_prop=p_prop, 
                min_seed=min_seed, 
                r_grow=r_grow, 
                tree_carrying_capacity=tree_carrying_capacity,
                neighborhood_carrying_capacity=n_cc, 
                max_ignite=max_ignite,
                rng=rng)

        
        ratio[trial_counter] = r
        #### BIOMASS
        # calculate summary stats for that trial and record them
        biomass = results_dict['biomass'][burnout_time:]
        # mean (total biomass / L*L)
        mean_biomass[trial_counter] = np.mean(biomass/(L*L))
        # median (total biomass / L*L)
        median_biomass[trial_counter] = np.median(biomass/(L*L))
        # std (total biomass / L*L)
        std_biomass[trial_counter] = np.std(biomass/(L*L))
        # 
        max_biomass[trial_counter] = np.max(biomass/(L*L))
        min_biomass[trial_counter] = np.min(biomass/(L*L))
        
        #### PATCH SIZE
        num_patches = results_dict['num_patches'][burnout_time:]
        avg_patch_size = results_dict['avg_patch_size'][burnout_time:]
        med_patch_size = results_dict['med_patch_size'][burnout_time:]
        
        # Get indices that are not nan so stuff doesn't break
        finite_idx = ~np.isnan(avg_patch_size)
        # mean (num patches)
        mean_patch_count[trial_counter] = np.nanmean(num_patches)
        # mean of the means
        mean_mean_patch_size[trial_counter] = np.nanmean(avg_patch_size)
        # sd of means
        std_mean_patch_size[trial_counter] = np.nanstd(avg_patch_size)
        # slope and int of mean patch size vs num patches

        # try to get slope, intercept
        # return nan if doesn't work 
        try: 
            slope, intercept = np.polyfit(np.log(num_patches[finite_idx]), np.log(avg_patch_size[finite_idx]), 1)
            ln_num_avg_patch_int[trial_counter] = intercept
            ln_num_avg_patch_slope[trial_counter] = slope
        except: 
            ln_num_avg_patch_int[trial_counter] = np.nan
            ln_num_avg_patch_slope[trial_counter] = np.nan

        trial_counter += 1
        

# Save results
df = pd.DataFrame({
    "ratio": ratio, 
    "mean_biomass": mean_biomass, 
    "median_biomass": median_biomass, 
    "max_biomass": max_biomass, 
    "min_biomass": min_biomass, 
    "std_biomass": std_biomass, 
    "mean_patch_count": mean_patch_count, 
    "mean_mean_patch_size": mean_mean_patch_size, 
    "std_mean_patch_size": std_mean_patch_size, 
    "ln_num_avg_patch_int": ln_num_avg_patch_int, 
    "ln_num_avg_patch_slope": ln_num_avg_patch_slope
})
df.to_csv("outputs/neighborcc_test.csv")


# Choose one reasonable value of neighborhood CC ratio for all other experiments
# Choose 5. 
neighborhood_carrying_capacity = tree_carrying_capacity*5

# growth rate
# run from 0.01 to 1 (intervals of 0.05)

ratios = range(20)
num_trials = len(ratios) * num_reps
ratio = np.zeros(num_trials)
mean_biomass = np.zeros(num_trials)
median_biomass = np.zeros(num_trials)
max_biomass = np.zeros(num_trials)
min_biomass = np.zeros(num_trials)
std_biomass = np.zeros(num_trials)
mean_patch_count = np.zeros(num_trials)
mean_mean_patch_size = np.zeros(num_trials)
std_mean_patch_size = np.zeros(num_trials)
ln_num_avg_patch_int = np.zeros(num_trials)
ln_num_avg_patch_slope = np.zeros(num_trials)

# keep track of trial number
trial_counter = 0

for r in ratios: 
    r = r + 1
    n_cc = tree_carrying_capacity*r
    for rep in range(num_reps):
        print(trial_counter)
        results_dict = run_simulation(m=m, 
                L=L, 
                t_steps=time_steps, 
                d=d, 
                init_grass=init_grass, 
                init_tree=init_tree, 
                p_disp=p_disp, 
                p_prop=p_prop, 
                min_seed=min_seed, 
                r_grow=r_grow, 
                tree_carrying_capacity=tree_carrying_capacity,
                neighborhood_carrying_capacity=n_cc, 
                max_ignite=max_ignite,
                rng=rng)

        
        ratio[trial_counter] = r
        #### BIOMASS
        # calculate summary stats for that trial and record them
        biomass = results_dict['biomass'][burnout_time:]
        # mean (total biomass / L*L)
        mean_biomass[trial_counter] = np.mean(biomass/(L*L))
        # median (total biomass / L*L)
        median_biomass[trial_counter] = np.median(biomass/(L*L))
        # std (total biomass / L*L)
        std_biomass[trial_counter] = np.std(biomass/(L*L))
        # 
        max_biomass[trial_counter] = np.max(biomass/(L*L))
        min_biomass[trial_counter] = np.min(biomass/(L*L))
        
        #### PATCH SIZE
        num_patches = results_dict['num_patches'][burnout_time:]
        avg_patch_size = results_dict['avg_patch_size'][burnout_time:]
        med_patch_size = results_dict['med_patch_size'][burnout_time:]
        
        # Get indices that are not nan so stuff doesn't break
        finite_idx = ~np.isnan(avg_patch_size)
        # mean (num patches)
        mean_patch_count[trial_counter] = np.nanmean(num_patches)
        # mean of the means
        mean_mean_patch_size[trial_counter] = np.nanmean(avg_patch_size)
        # sd of means
        std_mean_patch_size[trial_counter] = np.nanstd(avg_patch_size)
        # slope and int of mean patch size vs num patches

        # try to get slope, intercept
        # return nan if doesn't work 
        try: 
            slope, intercept = np.polyfit(np.log(num_patches[finite_idx]), np.log(avg_patch_size[finite_idx]), 1)
            ln_num_avg_patch_int[trial_counter] = intercept
            ln_num_avg_patch_slope[trial_counter] = slope
        except: 
            ln_num_avg_patch_int[trial_counter] = np.nan
            ln_num_avg_patch_slope[trial_counter] = np.nan

        trial_counter += 1


# carrying capacity
# run from 50 to 1000 (intervals of 50)
# leave min_seed consistent at 10
# leave growth rate consistent (choose a value from growth rate experiment)
# leave ratio of tree_carrying_capacity/neighborhood_carrying_capacity the same

# ignition probability
# run from 0.000001 to 0.1 (double for each iteration)