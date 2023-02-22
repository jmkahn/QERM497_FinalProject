'''this file is not super organized. contains snippets of code that we use 
to test the simulation''' 

# %%
from simulation import *
from viz import plot_slices_of_simulation
import matplotlib.pyplot as plt
output =  run_simulation(m=0.5, L=40, t_steps=5, d=1, init_grass=0.05, init_tree=0.95, p_ig_gmax=0.001, p_ig_tmax=0.0002, r_spr_tmax=0.8, r_spr_gmax=0.5, r_cat_tmax=0.3, r_cat_gmax=0.8, output_times=[0, 1, 2, 3, 4])
#plot_slices_of_simulation(output, [0,1,2,3,4])

#%%
from simulation import *
forest = initialize_forest(L=40, d=1, init_grass=0.05, init_tree=0.95)
params_dict = initialize_params_dict(m=0.5, L=40, t_steps=20, d=1, init_grass=0.05, init_tree=0.95, p_ig_gmax=0.003, p_ig_tmax=0.0001, r_spr_tmax=0.8, r_spr_gmax=.5, r_cat_tmax=0.3, r_cat_gmax=.8)
fire_season(forest, params_dict)

# %% 
# test set_probabilities function (w/out moisture dependence)
set_probabilities(0.5, 0.01, 0.002, 0.7, 0.5, 0.3, 0.7)
# return values should be 
# {"p_gro_ag": 0.5, "p_gro_gt": 0.5, "p_ig_g": 0.01, "p_ig_t": 0.002, "p_spr_gg": 0.35, "p_spr_tg": 0.49, "p_spr_gt": 0.15, "p_spr_tt": 0.21}

# %% 
# test run_simulation 
output = run_simulation(m=0.8, L=10, t_steps=5, d=2, init_grass=0.1, init_tree=0.1, p_ig_gmax=0.1, p_ig_tmax=0.02, r_spr_tmax=0.7, r_spr_gmax=0.5, r_cat_tmax=0.3, r_cat_gmax=0.7, output_times=[0, 1, 2, 3, 4])
plot_slices_of_simulation(output, [0,1,2,3,4])
# 
# %%
import matplotlib.colors as colors
colors_list = ["black", "yellow", "green"]
cmap = colors.ListedColormap(colors_list)
bounds = [0,1,2,3]
norm = colors.BoundaryNorm(bounds, cmap.N)

forest = initialize_forest(L=10, d=2, init_grass=0.2, init_tree=0.3)
plt.imshow(forest,  cmap=cmap, norm=norm)

params_dict = initialize_params_dict(m=0.5, L=10, t_steps=5, d=2, init_grass=0.3, init_tree=0.4, p_ig_gmax=0.01, p_ig_tmax=0.002, r_spr_tmax=0.7, r_spr_gmax=0.5, r_cat_tmax=0.3, r_cat_gmax=0.7)
forest = grow_season(forest, params_dict)
plt.imshow(forest,  cmap=cmap, norm=norm)

# %%
