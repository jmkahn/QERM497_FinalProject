from matplotlib import colors
import matplotlib.pyplot as plt
def plot_slices_of_simulation(slices, times, burn_masks=None):
    '''Params: 
        save_times: a list or np array of time steps 
        slices: an array with shape (2*save_times, _, _), where L is the side length of a 2-D world. 
            Each subarray contains the state of the forest pre- and post-fire at a time in save_times
        burn_mask (optional): array of binary masks, if included, burned areas will be displayed in orange
    Graphs the state of the simulation in a series of subplots, showing its progression over time 
    '''
    NUM_ROWS = 3
    colors_list = ["black", "yellow", "green", "orange"]
    cmap = colors.ListedColormap(colors_list)
    bounds = [0, 1, 2, 3, 4]
    norm = colors.BoundaryNorm(bounds, cmap.N)

    fig, ax = plt.subplots(nrows=NUM_ROWS, ncols=len(times), figsize=(20, 8))
    ax = ax.T.reshape(-1) # go down columns first
    for i in range(len(ax)):
        if i != 0: # only show labels for first subplot
            ax[i].xaxis.set_visible(False) 
            ax[i].yaxis.set_visible(False)
        ax[i].imshow(slices[i], cmap=cmap, norm=norm)

        # add burn_mask to during fire plots
        if (burn_masks is not None) and (i % NUM_ROWS == 1): 
            current_slice = slices[i] + 4*burn_masks[int(i/NUM_ROWS)]# 4 is the color that will be orange 
            ax[i].imshow(current_slice, cmap, norm)

        # and show 

        # label top row only 
        if (i % NUM_ROWS == 0): 
            ax[i].title.set_text("year = "+str(times[int(i/NUM_ROWS)]))

    # label top row as "pre-fire" and bottom row as "post-fire" #TODO: not done, looks bad rn 
    ax[0].set_ylabel("Pre-fire season", rotation = 0, fontsize = 20, labelpad = 50)


    fig.suptitle("State of simulation at various time steps")
    plt.show()
    
