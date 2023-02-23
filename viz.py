from matplotlib import colors
import matplotlib.pyplot as plt
def plot_slices_of_simulation(slices, times, burn_masks=None):
    '''Params: 
        times: a list of time steps 
        slices: an array with shape (2*save_times, _, _), where L is the side length of a 2-D world. 
            Each subarray contains the state of the forest pre- and post-fire at a time in save_times
        burn_mask (optional): array of binary masks, if included, burned areas will be displayed in orange
    Graphs the state of the simulation in a series of subplots, showing its progression over time 
    '''
    if burn_masks is not None: 
        num_rows = 3
    else: 
        num_rows = 2

    colors_list = ["black", "yellow", "green", "orange"]
    cmap = colors.ListedColormap(colors_list)
    bounds = [0, 1, 2, 3, 4]
    norm = colors.BoundaryNorm(bounds, cmap.N)

    fig, ax = plt.subplots(nrows=num_rows, ncols=len(times), figsize=(20, 8))
    for i in range(len(times)): 
        if i != 0: # only show labels for first subplot
            ax[0][i].xaxis.set_visible(False) 
            ax[0][i].yaxis.set_visible(False)
            ax[num_rows-1][i].xaxis.set_visible(False) 
            ax[num_rows-1][i].yaxis.set_visible(False)

        # show the plot after growth season 
        ax[0][i].imshow(slices[2*i], cmap=cmap, norm=norm)

        # add burn_mask to during fire plots
        if (burn_masks is not None): 
            with_mask = slices[2*i+1] + 4*burn_masks[i]# 4 is the color that will be orange 
            ax[1][i].imshow(with_mask, cmap, norm)
            ax[1][i].xaxis.set_visible(False) 
            ax[1][i].yaxis.set_visible(False)

        # and show post-fire without mask 
        ax[num_rows-1][i].imshow(slices[2*i+1], cmap, norm)


        # label top row only 
        ax[0][i].title.set_text("year = "+str(times[int(i/num_rows)]))

    # label top row as "pre-fire" and bottom row as "post-fire" #TODO: not done, looks bad rn 
    # ax[0].set_ylabel("Pre-fire season", rotation = 0, fontsize = 20, labelpad = 50)


    fig.suptitle("State of simulation at various time steps")
    plt.show()
    
