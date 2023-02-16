'''
This module gives the code to simulate a forest fire [add description]
Written by J Kahn and Helen Miller 
'''
import numpy as np 


def set_probabilities(m, p_ig_gmax, p_ig_tmax, r_spr_tmax, r_spr_gmax, r_cat_tmax, r_cat_gmax): #TODO: J
    '''
    Creates a dictionary of growth, spontaneous ignition, and fire spread probabilities for grass and trees based the input parameters and moisture levels
    Params: 
        - m : (float) Moisture level (fundamental parameter dictating probabilities of growth and fire spread) (can vary between 0 and 1, where 0 is no moisture and 1 is total saturation) 
        - p_ig_gmax : (float) Maximum probability of spontaneous grass ignition (when moisture=0) (p_ig_tmax < p_ig_gmax)
        - p_ig_tmax : (float) Maximum probability of spontaneous tree ignition (when moisture=0) (p_ig_tmax < p_ig_gmax)

    Generates the following parameters as a function of the input parameters: 
        - p_gro_ag: probability that grass propagates to neighboring ash site, equal to m 
        - p_gro_gt: probability that tree propagates to neighboring grass site, equal to m
        - p_ig_t: (prob of tree ignition) = p_ig_tmax*(1 - m)
        - p_ig_g: (prob of grass ignition) = p_ig_gmax*(1 - m)
        - p_spr_gg: prob of fire spreading from grass to neighboring grass
        - p_spr_tg: prob of fire spreading from tree to neighboring grass
        - p_spr_gt: prob of fire spreading from grass to neighboring tree
        - p_spr_tt: prob of fire spreading from tree to neighboring tree
    NOTE: probilities are named with the convention "p_[process]_[from state to state]"
    '''
    # Growth rates vary proportionately with moisture level  
    p_gro_ag = m
    p_gro_gt = m

    # Ignition probabilities vary inversely with moisture level 
    p_ig_g = p_ig_gmax #*(1 - m)
    p_ig_t = p_ig_tmax #*(1 - m) #TODO: add in moisture dependence later 

    # Fire spread between neighbors depends on moisture and state each cell 
    # Multiply (rate of spread) times (rate of catching) for each combination of (grass, tree) to get probabilities of fire spreading between trees and grass 
    p_spr_gg = r_spr_gmax* r_cat_gmax #*(1-m)**2 #TODO: add moisture dependence 
    p_spr_tg = r_spr_tmax* r_cat_gmax #*(1-m)**2
    p_spr_gt = r_spr_gmax* r_cat_tmax #*(1-m)**2
    p_spr_tt = r_spr_tmax* r_cat_tmax #*(1-m)**2

    # Returns an dict of all the transition probabilities
    return {"p_gro_ag": p_gro_ag, "p_gro_gt": p_gro_gt, "p_ig_g": p_ig_g, "p_ig_t": p_ig_t, "p_spr_gg": p_spr_gg, "p_spr_tg": p_spr_tg, "p_spr_gt": p_spr_gt, "p_spr_tt": p_spr_tt}

    # %% 
    # test set_probabilities function (w/out moisture dependence)
    set_probabilities(0.5, 0.01, 0.002, 0.7, 0.5, 0.3, 0.7)
    # return values should be 
    # {"p_gro_ag": 0.5, "p_gro_gt": 0.5, "p_ig_g": 0.01, "p_ig_t": 0.002, "p_spr_gg": 0.35, "p_spr_tg": 0.49, "p_spr_gt": 0.15, "p_spr_tt": 0.21}


def initialize_forest(L, d, init_grass, init_tree): 
    ''' initialize_forest creates the world for the simulation (stored in a 2D array) and populates it with a randomly dispersed initial set of trees and grass. 
    Ash, grass and trees are represented by ints: 0=ash; 1=grass; 2=tree.'''
    # create (L+2d)x(L+2d) array(forest)
    # Randomly populate forest with trees and grass based on initial conditions
    forest = np.random.choice([0,1,2],
                              size=(L+2*d,L+2*d), 
                              p=[1-(init_tree+init_grass), init_grass, init_tree])
    # 0=ash; 1=grass; 2=tree
    # Set boundary region (2d length around edge) to be all ash 
    forest[:, 0:d] = 0
    forest[:, (L+d):] = 0
    forest[0:d, :] = 0
    forest[(L+d):, :] = 0
    
    return forest


def grow_season(world, probs_dict): #TODO: J 
    '''grow_season runs one iteration of the world in which new vegetation grows'''
    # Iterate across the forest grid 
    # For each cell: 
    # Get neighborhood 
    # Count number of each state in neighborhood
    # Get its current state 
    # If ash: roll to see if grows grass or stays ash
    # If grass: roll to see if grows tree or stays grass 
    # ^(these depend on number of neighbors in grass/tree state)
    # Return world at next time step (and optional data) 
    return world


def fire_season(world, probs_dict): #TODO: H
    '''fire_season runs one iteration of the world, in which wildfires initiate, propagate, and die out'''
    # Create two lists: 
    # Fire_location = list of locations of currently on-fire cells in grid
    # Fire_type = “T” or “G” depending on whether each cell is tree or grass fire (aligned with Fire_location list) 
    # Randomly select some of the vegetation cells to ignite and add to the lists 
    # While there are still cells on fire:
    # Take the first element in fire_location: 	
    # Get its neighbors (3x3 moore neighborhood)
    # Randomly select which of its neighbors get set on fire (probabilities depend on moisture and focal/neighbor state)
    # Add location of each neighbor which gets set on fire gets to end of fire_location list
    # Add type of each neighbor which gets set on fire to fire_type end of list
    # Update world
    # Change focal cell to ash
    # Change appropriate neighbors to fire
    # Remove first element from fire_location list
    # Remove first element from fire_type list
    # (Optionally record data)
    # Return world at next time step (and optional data) 
    return world

def run_simulation(m, L, t_steps, d, init_grass, init_tree, p_ig_gmax, p_ig_tmax, r_spr_tmax, r_spr_gmax, r_cat_tmax, r_cat_gmax, output_times=[]): 
    ''' run_simulation is the wrapper function which initializes the simulation and runs it for a specified number of growth and fire seasons
    Params:  
        m : (float) Moisture level (fundamental parameter dictating probabilities of growth and fire spread) (can vary between 0 and 1, where 0 is no moisture and 1 is total saturation) 
        L : (int) side length of forest
        t_steps : (int) Number of “years” (alternations between growth and fire season) to run simulation for
        d : (int) Radius of neighborhood during growth season 
        p_ig_gmax : (float) Maximum probability of spontaneous grass ignition (when moisture=0) (p_ig_tmax < p_ig_gmax)
        p_ig_tmax : (float) Maximum probability of spontaneous tree ignition (when moisture=0) (p_ig_tmax < p_ig_gmax)
        r_spr_tmax : (float) Max rate at which fire spreads from a tree (r_spr_gmax < r_spr_tmax) 
        r_spr_gmax : (float) Max rate at which fire spreads from grass to neighbors (r_spr_gmax < r_spr_tmax
        r_cat_tmax : (float) Max rate at which trees catch fire from nearby source (r_cat_tmax < r_cat_gmax)
        r_cat_gmax : (float) Max rate at which grass catches fire from nearby source (r_cat_tmax < r_cat_gmax)
        init_tree : (float) Initial proportion of tree cover (0 to 1; init_tree + init_grass <= 1) 
        init_grass : (float) Initial proportion of grass cover (0 to 1; init_tree + init_grass <= 1)
    '''
    forest = initialize_forest(L, d, init_grass, init_tree)
    probs_dict = set_probabilities(m, p_ig_gmax, p_ig_tmax, r_spr_tmax, r_spr_gmax, r_cat_tmax, r_cat_gmax) # Set all the probabilities based on moisture levels 

    save_counter = 0 # save state data at select times 
    output_slices = np.zeros((len(output_times), (L+2*d), (L+2*d)))
    for t in range(t_steps): 
        forest = grow_season(forest, probs_dict) #TODO: return data 
        forest = fire_season(forest, probs_dict)
        if t in output_times: 
            output_slices[save_counter] = forest
            save_counter += 1
    return output_slices

# %% 
# test run_simulation 
run_simulation(m=0.5, L=10, t_steps=5, d=5, init_grass=0.3, init_tree=0.4, p_ig_gmax=0.01, p_ig_tmax=0.002, r_spr_tmax=0.7, r_spr_gmax=0.5, r_cat_tmax=0.3, r_cat_gmax=0.7, output_times=[1, 2])
# 