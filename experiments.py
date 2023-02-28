"""Run experiments over a range of parameters of interest"""

# TODO: run some parameter sweeps on: 
    # * neighborhood CC (in relationship to tree cc)
    # * growth rate
    # * carrying capacity
    # * ignition probability
    # collect data on mean and sd of biomass and patch size (from 200-100 years or something)
import numpy as np
from simulation import run_simulation

# DEFAULT PARAMS
m = 0.5 # placeholder, not used. 
L = 100
time_steps = 500
d = 4
init_grass=0.2 
init_tree=0.1
p_disp=.05
p_prop=0.1
min_seed=10
r_grow=0.2
tree_carrying_capacity = 100
neighborhood_carrying_capacity = 500
max_ignite=0.01
rand_seed = 12345

num_trials = 10

# Neighborhood CC 
# Run over a range of different ratios of tree_carrying_capacity/neighborhood_carrying_capacity
# from 1/10 to 1 (1/10, 1/9, 1/8, etc)

# set up data collection
# will have a couple of summary stats for each iteration


for i in range(10): 
    i = i + 1
    n_cc = tree_carrying_capacity*i
    for trial in range(num_trials):
        r_seed = rand_seed + trial
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
                rand_seed=r_seed)
        # calculate summary stats for that trial and record them

# Save results



# Choose one reasonable value of neighborhood CC ratio for all other experiments

# growth rate
# run from 0.01 to 1 (intervals of 0.5)


# carrying capacity
# run from 50 to 1000 (intervals of 50)
# leave min_seed consistent at 10
# leave growth rate consistent (choose a value from growth rate experiment)
# leave ratio of tree_carrying_capacity/neighborhood_carrying_capacity the same

# ignition probability
# run from 0.000001 to 0.1 (double for each iteration)