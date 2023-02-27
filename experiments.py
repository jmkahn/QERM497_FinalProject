# Note: this can be exported as jupyter notebook by
# Ctl+Shift+P > Jupyter: Export Current Python File as Jupyter Notebook

# %%
from simulation import *
from viz import plot_slices_of_simulation, veg_fire_over_time, plot_slices_of_simulation_biomass
import matplotlib.pyplot as plt

# %% 
# Test with new biomass version
# Params
m = 0.5 # placeholder, not used. 
L = 10
time_steps = 10
d = 2
init_grass=0.4 
init_tree=0.5
p_disp=.01
p_prop=0.1
min_seed=20
r_grow=0.01
tree_carrying_capacity = 10
neighborhood_carrying_capacity = 150
max_ignite=0.1
output_times=[0,2,4,6,8]

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
               neighborhood_carrying_capacity=neighborhood_carrying_capacity, 
               max_ignite=max_ignite, 
               output_times=output_times)
plot_slices_of_simulation_biomass(results_dict=results_dict, plot_burn_masks=False)

# %%
slices = results_dict['output_slices']


# %%
# results_dict = run_simulation(m=0.4, L=50, t_steps=50, d=2, init_grass=0.1, init_tree=0.1, p_gro_gmax = 0.02, p_lightning=0.01, r_spr_tmax=1, r_spr_gmax=0.8, r_cat_tmax=0.4, r_cat_gmax=1, output_times=[0, 10, 20, 30, 49])
# plot_slices_of_simulation(results_dict["output_slices"] , [0,1,2,3,4], burn_masks=results_dict["burned_area_masks"])

# results_dict = run_simulation(m=0.2, L=50, t_steps=50, d=2, init_grass=0.1, init_tree=0.1, p_gro_gmax = 0.02, p_lightning=0.01, r_spr_tmax=1, r_spr_gmax=0.8, r_cat_tmax=0.4, r_cat_gmax=1, output_times=[0, 10, 20, 30, 49])
# plot_slices_of_simulation(results_dict["output_slices"] , [0,1,2,3,4], burn_masks=results_dict["burned_area_masks"])

# results_dict = run_simulation(m=0.1, L=50, t_steps=50, d=2, init_grass=0.1, init_tree=0.1, p_gro_gmax = 0.02, p_lightning=0.01, r_spr_tmax=1, r_spr_gmax=0.8, r_cat_tmax=0.4, r_cat_gmax=1, output_times=[0, 10, 20, 30, 49])
# plot_slices_of_simulation(results_dict["output_slices"] , [0,1,2,3,4], burn_masks=results_dict["burned_area_masks"])

# %% [Markdown]


# %%
# Initialize parameters
p_disp = 0.01   # p Grass grows in random ash cell (long-range dispersal)
# OR random grass cell becomes tree
p_gro_ag = 0.3  # p grass propogates to neighboring ash cell
# p tree propogates to neighboring cell (neighboring ash cells become grass)
p_gro_gt = 0.2
p_ig_g = 0.001  # p any given grass cell ignites randomly at beginning of season
p_ig_t = 0.0005  # p any given tree cell ignites randomly at beginning of season
p_spr_gg = 0.25  # p grass fire spreads to neighbor grass
p_spr_tg = 0.25  # p tree fire spreads to neighbor grass
p_spr_gt = 0.25  # p grass fire spreads to neighbor tree
p_spr_tt = 0.25  # p tree fire spreads to neighbor tree

L = 50
time_steps = 100
output_times = [0,1,2,30,49]

# %%
results_dict = run_simulation(m=0.5, # not currently used
                              L=L,
                              t_steps=time_steps,
                              d=2,
                              init_grass=0.1,
                              init_tree=0.1,
                              p_disp=p_disp,
                              p_gro_ag=p_gro_ag,
                              p_gro_gt=p_gro_gt,
                              p_ig_g=p_ig_g,
                              p_ig_t=p_ig_t,
                              p_spr_gg=p_spr_gg,
                              p_spr_tg=p_spr_tg,
                              p_spr_gt=p_spr_gt,
                              p_spr_tt=p_spr_tt,
                              output_times=output_times)

fig, ax = plot_slices_of_simulation(results_dict)
fig.savefig("plots/defaults_maps.png")
fig.show()

fig = veg_fire_over_time(time_steps, results_dict)
fig.savefig("plots/defaults_timeseries.png")



# %% [markdown] 
# Now the only thing to do, I think, is to start changing variables one at a time 
# and save the outputs, and see if we can come up with some interesting dynamics, 
# and get an idea of how parameters interact. 

# %%
# 
# For grass-grass fire:  
for r in np.linspace(0, 1, 11): 
    results_dict = run_simulation(m=0.5, # not currently used
                              L=50,
                              t_steps=time_steps,
                              d=2,
                              init_grass=0.1,
                              init_tree=0.1,
                              p_disp=p_disp,
                              p_gro_ag=p_gro_ag,
                              p_gro_gt=p_gro_gt,
                              p_ig_g=p_ig_g,
                              p_ig_t=p_ig_t,
                              p_spr_gg=r, # TESTING THIS ONE
                              p_spr_tg=p_spr_tg,
                              p_spr_gt=p_spr_gt,
                              p_spr_tt=p_spr_tt,
                              output_times=output_times)
    r_txt = str(round(r, 5))
    fig, ax = plot_slices_of_simulation(results_dict, title="p_spr_gg="+r_txt)
    fig.savefig("plots/p_spr_gg-"+r_txt+"-maps.png")
    fig = veg_fire_over_time(time_steps, results_dict, title="p_spr_gg="+r_txt)
    fig.savefig("plots/p_spr_gg-"+r_txt+"-timeseries.png")

# tree - tree fire
for r in np.linspace(0, 1, 11): 
    results_dict = run_simulation(m=0.5, # not currently used
                              L=50,
                              t_steps=time_steps,
                              d=2,
                              init_grass=0.1,
                              init_tree=0.1,
                              p_disp=p_disp,
                              p_gro_ag=p_gro_ag,
                              p_gro_gt=p_gro_gt,
                              p_ig_g=p_ig_g,
                              p_ig_t=p_ig_t,
                              p_spr_gg=p_spr_gg,
                              p_spr_tg=p_spr_tg,
                              p_spr_gt=p_spr_gt,
                              p_spr_tt=r, # TESTING THIS ONE
                              output_times=output_times)
    r_txt = str(round(r, 5))
    fig, ax = plot_slices_of_simulation(results_dict, title="p_spr_tt="+r_txt)
    fig.savefig("plots/p_spr_tt-"+r_txt+"-maps.png")
    fig = veg_fire_over_time(time_steps, results_dict, title="p_spr_tt="+r_txt)
    fig.savefig("plots/p_spr_tt-"+r_txt+"-timeseries.png")


# tree - grass fire
for r in np.linspace(0, 1, 11): 
    results_dict = run_simulation(m=0.5, # not currently used
                              L=50,
                              t_steps=time_steps,
                              d=2,
                              init_grass=0.1,
                              init_tree=0.1,
                              p_disp=p_disp,
                              p_gro_ag=p_gro_ag,
                              p_gro_gt=p_gro_gt,
                              p_ig_g=p_ig_g,
                              p_ig_t=p_ig_t,
                              p_spr_gg=p_spr_gg,
                              p_spr_tg=r, # TESTING THIS ONE
                              p_spr_gt=p_spr_gt,
                              p_spr_tt=p_spr_tt,
                              output_times=output_times)
    r_txt = str(round(r, 5))
    fig, ax = plot_slices_of_simulation(results_dict, title="p_spr_tg="+r_txt)
    fig.savefig("plots/p_spr_tg-"+r_txt+"-maps.png")
    fig = veg_fire_over_time(time_steps, results_dict, title="p_spr_tg="+r_txt)
    fig.savefig("plots/p_spr_tg-"+r_txt+"-timeseries.png")

# grass - tree fire 
for r in np.linspace(0, 1, 11): 
    results_dict = run_simulation(m=0.5, # not currently used
                              L=50,
                              t_steps=time_steps,
                              d=2,
                              init_grass=0.1,
                              init_tree=0.1,
                              p_disp=p_disp,
                              p_gro_ag=p_gro_ag,
                              p_gro_gt=p_gro_gt,
                              p_ig_g=p_ig_g,
                              p_ig_t=p_ig_t,
                              p_spr_gg=p_spr_gg,
                              p_spr_tg=p_spr_tg,
                              p_spr_gt=r,# TESTING THIS ONE
                              p_spr_tt=p_spr_tt,
                              output_times=output_times)
    r_txt = str(round(r, 5))
    fig, ax = plot_slices_of_simulation(results_dict, title="p_spr_gt="+r_txt)
    fig.savefig("plots/p_spr_gt-"+r_txt+"-maps.png")
    fig = veg_fire_over_time(time_steps, results_dict, title="p_spr_gt="+r_txt)
    fig.savefig("plots/p_spr_gt-"+r_txt+"-timeseries.png")