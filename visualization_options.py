def plot_slices_of_simulation(slices, times, neighborhood_size):
    '''Params: 
        save_times: a list or np array of time steps 
        slices: an array with shape (save_times, L, L), where L is the side length of a 2-D world. Each subarray contains the state of the world at a time in save_times

    Graphs the state of the simulation in a series of subplots, showing its progression over time 
    '''
    fig, ax = plt.subplots(nrows=1, ncols=len(save_times), figsize=(20, 4))
    ax = ax.reshape(-1)
    for i in range(len(ax)):
        ax[i].xaxis.set_visible(False) 
        ax[i].yaxis.set_visible(False)
        ax[i].imshow(slices[i])
        ax[i].title.set_text("t="+str(save_times[i]))
    fig.suptitle("State of simulation at various time steps")
    plt.show()
