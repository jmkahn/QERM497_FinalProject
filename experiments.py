# Note: this can be exported as jupyter notebook by
# Ctl+Shift+P > Jupyter: Export Current Python File as Jupyter Notebook

# %%
from simulation import *
from viz import plot_slices_of_simulation, veg_fire_over_time, plot_slices_of_simulation_biomass
import matplotlib.pyplot as plt
from matplotlib import colors

def plot_slices_of_simulation_biomass(results_dict, title = "", plot_burn_masks=True): 
    '''Params: 
    times: a list of time steps 
    slices: an array with shape (2*save_times, _, _), where L is the side length of a 2-D world. 
        Each subarray contains the state of the forest pre- and post-fire at a time in save_times
    burn_mask (optional): array of binary masks, if included, burned areas will be displayed in orange
    Graphs the state of the simulation in a series of subplots, showing its progression over time 
    '''
    slices = results_dict["output_slices"]
    times = results_dict["output_times"]
    tree_cc = results_dict["params_dict"]["tree_carrying_capacity"]
    if plot_burn_masks: 
        num_rows = 3
        burn_masks = results_dict["burned_area_masks"]
    else: 
        num_rows = 2

    fig, ax = plt.subplots(nrows=num_rows, ncols=len(times), figsize=(20, 8))
    for i in range(len(times)): 
        ax[0][i].xaxis.set_visible(False) #hide axis ticks 
        ax[0][i].yaxis.set_visible(False)
        ax[num_rows-1][i].xaxis.set_visible(False) 
        ax[num_rows-1][i].yaxis.set_visible(False)

        # show the plot after growth season 
        ax[0][i].imshow(slices[2*i], vmin=1, vmax=tree_cc*1.2)
        
        # Make all ash cells black
        colors_list = ["black"]
        cmap = colors.ListedColormap(colors_list)
        bounds = [0,1]
        norm = colors.BoundaryNorm(bounds, cmap.N)
        forest_ash_before = (slices[2*i] == 0).astype(float)
        ax[0][i].imshow(forest_ash_before, alpha=forest_ash_before, cmap = cmap, norm = norm)

        
        # add burn_mask to during fire plots
        if (plot_burn_masks): 
            
            ax[1][i].imshow(slices[2*i+1],  vmin=1, vmax=tree_cc*1.2)
            ax[1][i].imshow(forest_ash_before, alpha=forest_ash_before, cmap = cmap, norm = norm)
            cmap_burn = colors.ListedColormap(["black", "orange"])
            bounds_burn = [0,1,2]
            norm_burn = colors.BoundaryNorm(bounds_burn, cmap_burn.N)
            ax[1][i].imshow(burn_masks[i], 
                            alpha=burn_masks[i].astype(float), 
                            cmap=cmap_burn, 
                            norm=norm_burn)
            ax[1][i].xaxis.set_visible(False) 
            ax[1][i].yaxis.set_visible(False)

        # and show post-fire without mask 
        ax[num_rows-1][i].imshow(slices[2*i], vmin=1, vmax=tree_cc*1.2,)
        forest_ash_after = (slices[2*i+1] == 0).astype(float)
        ax[num_rows-1][i].imshow(forest_ash_after, alpha=forest_ash_after, cmap = cmap, norm = norm)


        # label top row only 
        ax[0][i].title.set_text("year = "+str(times[i]))

    # label top row as "pre-fire" and bottom row as "post-fire" #TODO: not done, looks bad rn 
    # ax[0].set_ylabel("Pre-fire season", rotation = 0, fontsize = 20, labelpad = 50)

    #only show axis ticks for first plot
    ax[0][0].xaxis.set_visible(True) 
    ax[0][0].yaxis.set_visible(True)
    if len(title) > 0: 
        plt.suptitle(title) 
    else: 
        fig.suptitle("State of simulation at various time steps")
    return fig, ax

from scipy import ndimage
def get_patch_sizes(bin_array, neighboorhood="moore"): 
    """Get an array of patch sizes from a binary array, where a patch 
    is a set of neighboring True cells. 
    Parameters: 
        - bin_array: binary array
        - neighborhood: neighborhood type. "moore" or "vneum" accepted"""
    
    # create structure object to know who to consider neighbors
    if neighboorhood == "moore": 
        s=ndimage.generate_binary_structure(2,2)
    elif neighboorhood == "vneum": 
        s=ndimage.generate_binary_structure(2,1)
    # Get an array where cells in each patch are represented by a different number
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.ndimage.label.html
    veg_patches, num_veg_patches = ndimage.label(bin_array, structure=s)
    # Get an array of length num_veg_patches where each element is the size 
    # of a patch
    patch_sizes = [np.sum(veg_patches == i) for i in range(1, num_veg_patches+1)]
    return patch_sizes

def plot_results(results_dict, title, filename):
    fig, ax = plot_slices_of_simulation_biomass(results_dict=results_dict, plot_burn_masks=True)
    fig.savefig("plots/"+filename+"-slices.png")

    fig = plt.figure()
    plt.scatter(range(time_steps), results_dict['area_burned'], color ="orange")
    plt.title(title + ": area burned per year")
    #plt.plot(range(time_steps), results_dict['biomass'], color = "green")
    fig.savefig("plots/"+filename+"-areaburned.png")

    fig = plt.figure()
    plt.plot(range(time_steps), results_dict['biomass'], color = "green")
    plt.title(title + ": total biomass over time")
    fig.savefig("plots/"+filename+"-biomass.png")


    
    slices = results_dict["output_slices"]
    times = results_dict["output_times"]
    fig, ax = plt.subplots(nrows=1, ncols=len(times), figsize=(25,4))
    for i, t in enumerate(times): 
        patch_sizes = get_patch_sizes(slices[2*i+1] > 0)
        patch_size_counts = np.unique(patch_sizes, return_counts=True)
        ax[i].scatter(np.log(patch_size_counts[0]), np.log(patch_size_counts[1]))
        ax[i].title.set_text("year = "+str(times[i]))
    fig.suptitle(title+": Distribution of patch sizes (log-log)")
    fig.savefig("plots/"+filename+"-patches.png")
# %% 
# Test with new biomass version
# Params
m = 0.5 # placeholder, not used. 
L = 100
time_steps = 500
d = 4
init_grass=0.2 
init_tree=0.1
p_disp=.05
p_prop=0.1
min_seed=20
r_grow=0.2
tree_carrying_capacity = 20
neighborhood_carrying_capacity = 90
max_ignite=0.1
output_times=[0,100,200,300,400,499]


results_dict_1 = run_simulation(m=m, 
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

plot_results(results_dict_1, title="experiment 1", filename="exp1")


# %%
# Same thing but with different params: larger carrying capacity

init_grass=0.2 
init_tree=0.1
p_disp=.05
p_prop=0.1
min_seed=20
r_grow=0.2
tree_carrying_capacity = 100
neighborhood_carrying_capacity = 500
max_ignite=0.1
output_times=[0,100,200,300,400,499]


results_dict_2 = run_simulation(m=m, 
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
plot_results(results_dict_2, title = "experiment 2", filename = "exp2")

# %%
# Again! But this time decrease ignition probability
m = 0.5 # placeholder, not used. 
L = 100
time_steps = 500
d = 4
init_grass=0.2 
init_tree=0.1
p_disp=.05
p_prop=0.1
min_seed=20
r_grow=0.2
tree_carrying_capacity = 20
neighborhood_carrying_capacity = 90
max_ignite=0.01
output_times=[0,100,200,300,400,499]


results_dict_3 = run_simulation(m=m, 
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
plot_results(results_dict_3, title="experiment 3", filename="exp3")

np.unique(results_dict_3['output_slices'][11][d:-d,d:-d].flatten(), return_counts=True)

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