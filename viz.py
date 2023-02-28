import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt
from scipy import ndimage

def plot_slices_of_simulation(results_dict, title = "", plot_burn_masks=True):
    '''Params: 
        times: a list of time steps 
        slices: an array with shape (2*save_times, _, _), where L is the side length of a 2-D world. 
            Each subarray contains the state of the forest pre- and post-fire at a time in save_times
        burn_mask (optional): array of binary masks, if included, burned areas will be displayed in orange
    Graphs the state of the simulation in a series of subplots, showing its progression over time 
    '''
    slices = results_dict["output_slices"]
    times = results_dict["output_times"]
    if plot_burn_masks: 
        num_rows = 3
        burn_masks = results_dict["burned_area_masks"]
    else: 
        num_rows = 2
    

    colors_list = ["black", "yellow", "green", "orange"]
    cmap = colors.ListedColormap(colors_list)
    bounds = [0, 1, 2, 3, 4]
    norm = colors.BoundaryNorm(bounds, cmap.N)

    fig, ax = plt.subplots(nrows=num_rows, ncols=len(times), figsize=(20, 8))
    for i in range(len(times)): 
        ax[0][i].xaxis.set_visible(False) #hide axis ticks 
        ax[0][i].yaxis.set_visible(False)
        ax[num_rows-1][i].xaxis.set_visible(False) 
        ax[num_rows-1][i].yaxis.set_visible(False)

        # show the plot after growth season 
        ax[0][i].imshow(slices[2*i], cmap=cmap, norm=norm)

        # add burn_mask to during fire plots
        if (plot_burn_masks): 
            with_mask = slices[2*i+1] + 4*burn_masks[i]# 4 is the color that will be orange 
            ax[1][i].imshow(with_mask, cmap, norm)
            ax[1][i].xaxis.set_visible(False) 
            ax[1][i].yaxis.set_visible(False)

        # and show post-fire without mask 
        ax[num_rows-1][i].imshow(slices[2*i+1], cmap, norm)


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

def veg_fire_over_time(t_steps, results_dict, title = ""):
    '''line plot of amount of grass, trees, and fire over time'''
    fig = plt.figure(figsize=(12, 6))
    plt.plot(range(t_steps), results_dict['grass_count'], color="yellow")
    plt.plot(range(t_steps), results_dict['tree_count'], color='green')
    plt.plot(range(t_steps), results_dict['area_burned'], color='orange')
    if len(title) > 0: 
        plt.suptitle(title) 
    else: 
        plt.suptitle("Dynamics of simulation over time")
    plt.legend(fontsize="large")
    plt.xlabel("Time steps ('years')")
    plt.ylabel("Units of Vegetation/Fire")
    return fig

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
            #ax[1][i].imshow(forest_ash_before, alpha=forest_ash_before, cmap = cmap, norm = norm)
            cmap_burn = colors.ListedColormap(["black", "orange"])
            bounds_burn = [0,1,2]
            norm_burn = colors.BoundaryNorm(bounds_burn, cmap_burn.N)
            ax[1][i].imshow(burn_masks[times[i]], 
                            alpha=burn_masks[times[i]].astype(float), 
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
