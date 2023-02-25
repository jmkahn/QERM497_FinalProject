'''this file is not super organized. contains snippets of code that we use 
to test the simulation'''

# %%
import matplotlib.colors as colors
from simulation import *
from viz import plot_slices_of_simulation
import matplotlib.pyplot as plt

# %%
output = run_simulation(m=0.5, L=40, t_steps=5, d=1, init_grass=0.05, init_tree=0.95, p_ig_gmax=0.001,
                        p_ig_tmax=0.0002, r_spr_tmax=0.8, r_spr_gmax=0.5, r_cat_tmax=0.3, r_cat_gmax=0.8, output_times=[0, 1, 2, 3, 4])
#plot_slices_of_simulation(output, [0,1,2,3,4])

# %%
# Initialize all parameters:
p_disp = 0.01   # p Grass grows in random ash cell (long-range dispersal)
# OR random grass cell becomes tree
p_gro_ag = 0.3  # p grass propogates to neighboring ash cell
# p tree propogates to neighboring cell (neighboring ash cells become grass)
p_gro_gt = 0.3
p_ig_g = 0.001  # p any given grass cell ignites randomly at beginning of season
p_ig_t = 0.001  # p any given tree cell ignites randomly at beginning of season
p_spr_gg = 0.25  # p grass fire spreads to neighbor grass
p_spr_tg = 0.25  # p tree fire spreads to neighbor grass
p_spr_gt = 0.25  # p grass fire spreads to neighbor tree
p_spr_tt = 0.25  # p tree fire spreads to neighbor tree
forest = initialize_forest(L=40, d=1, init_grass=0.05, init_tree=0.95)
params_dict = initialize_params_dict(m=0.5, L=40, t_steps=20, d=1, init_grass=0.05, init_tree=0.95,
                                     p_disp=p_disp,
                                     p_gro_ag=p_gro_ag,
                                     p_gro_gt=p_gro_gt,
                                     p_ig_g=p_ig_g,
                                     p_ig_t=p_ig_t,
                                     p_spr_gg=p_spr_gg,
                                     p_spr_tg=p_spr_tg,
                                     p_spr_gt=p_spr_gt,
                                     p_spr_tt=p_spr_tt)
forest, area_burned, burned_cells = fire_season(forest, params_dict)
plt.imshow(forest)
# %%
# test set_probabilities function (w/out moisture dependence)
set_probabilities(0.5, 0.01, 0.002, 0.7, 0.5, 0.3, 0.7)
# return values should be
# {"p_gro_ag": 0.5, "p_gro_gt": 0.5, "p_ig_g": 0.01, "p_ig_t": 0.002, "p_spr_gg": 0.35, "p_spr_tg": 0.49, "p_spr_gt": 0.15, "p_spr_tt": 0.21}

# %%
# test run_simulation
results_dict = run_simulation(m=0.4, L=50, t_steps=50, d=2, init_grass=0.1, init_tree=0.1, p_gro_gmax=0.02,
                              p_lightning=0.01, r_spr_tmax=1, r_spr_gmax=0.8, r_cat_tmax=0.4, r_cat_gmax=1, output_times=[0, 10, 20, 30, 49])
slices = plot_slices_of_simulation(results_dict["output_slices"], [
                          0, 1, 2, 3, 4], burn_masks=results_dict["burned_area_masks"])

# %%
# test run_simulation
t_steps = 10000
results_dict = run_simulation(m=0.6, L=100, t_steps=t_steps, d=5, init_grass=0.1, init_tree=0.1, p_ig_gmax=0.1, p_ig_tmax=0.02,
                              r_spr_tmax=0.03, r_spr_gmax=0.05, r_cat_tmax=0.03, r_cat_gmax=0.07, output_times=[0, 10, 20, 30, 40, 50, 60])
plot_slices_of_simulation(results_dict["output_slices"], [
                          0, 10, 20, 30, 40, 50, 60], burn_masks=results_dict["burned_area_masks"])

# %%
plt.figure(figsize=(12, 6))
plt.plot(range(t_steps), results_dict['grass_count'],
         color="yellow", label="Early Vegetation")
plt.plot(range(t_steps), results_dict['tree_count'],
         color='green', label="Mature Vegetation")
plt.plot(range(t_steps),
         results_dict['area_burned'], color='orange', label="Fire")
plt.suptitle("Dynamics of simulation over time")
plt.legend(fontsize="large")
plt.xlabel("Time steps ('years')")
plt.ylabel("Units of Vegetation/Fire")
# %%

colors_list = ["black", "yellow", "green"]
cmap = colors.ListedColormap(colors_list)
bounds = [0, 1, 2, 3]
norm = colors.BoundaryNorm(bounds, cmap.N)

forest = initialize_forest(L=10, d=2, init_grass=0.2, init_tree=0.3)
plt.imshow(forest,  cmap=cmap, norm=norm)

params_dict = initialize_params_dict(m=0.5, L=10, t_steps=5, d=2, init_grass=0.3, init_tree=0.4,
                                     p_ig_gmax=0.01, p_ig_tmax=0.002, r_spr_tmax=0.7, r_spr_gmax=0.5, r_cat_tmax=0.3, r_cat_gmax=0.7)
forest = grow_season(forest, params_dict)
plt.imshow(forest,  cmap=cmap, norm=norm)

# %%
t_steps = 100
m = 0.6
L = 100
d = 6
init_grass = 0.5
init_tree = 0.5
p_ig_gmax = 0.02
p_ig_tmax = 0.03
results_dict = run_simulation(m=m, L=L, t_steps=t_steps, d=d, init_grass=init_grass, init_tree=init_tree, p_ig_gmax=p_ig_gmax,
                              p_ig_tmax=p_ig_tmax, r_spr_tmax=0.03, r_spr_gmax=0.05, r_cat_tmax=0.03, r_cat_gmax=0.07, output_times=[0, 10, 20, 30, 40, 50, 60])
