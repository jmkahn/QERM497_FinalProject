# %%
from simulation import run_simulation
from visualization_options import plot_slices_of_simulation
output =  run_simulation(m=0.5, L=40, t_steps=5, d=1, init_grass=0.05, init_tree=0.95, p_ig_gmax=0.001, p_ig_tmax=0.0002, r_spr_tmax=0.8, r_spr_gmax=0.5, r_cat_tmax=0.3, r_cat_gmax=0.8, output_times=[0, 1, 2, 3, 4])
plot_slices_of_simulation(output, [0,1,2,3,4])
