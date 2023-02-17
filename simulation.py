'''
This module gives the code to simulate a forest fire [add description]
Written by J Kahn and Helen Miller 

Model Parameters:  
    L : (int) side length of forest
    t_steps : (int) Number of “years” (alternations between growth and fire season) to run simulation for
    d : (int) Radius of neighborhood during growth season 
    m : (float) Moisture level (fundamental parameter dictating probabilities of growth and fire spread) (can vary between 0 and 1, where 0 is no moisture and 1 is total saturation) 
    p_i_gmax : (float) Maximum probability of spontaneous grass ignition (when moisture=0) (p_i_tmax < p_i_gmax)
    p_i_tmax : (float) Maximum probability of spontaneous tree ignition (when moisture=0) (p_i_tmax < p_i_gmax)
    r_t_spread_max : (float) Max rate at which fire spreads from a tree (r_g_spread_max < r_t_spread_max) 
    r_g_spread_max : (float) Max rate at which fire spreads from grass to neighbors (r_g_spread_max < r_t_spread_max
    r_t_catch_max : (float) Max rate at which trees catch fire from nearby source (r_t_catch_max < r_g_catch_max)
    r_g_catch_max : (float) Max rate at which grass catches fire from nearby source (r_t_catch_max < r_g_catch_max)
    init_tree : (float) Initial proportion of tree cover (0 to 1; init_tree + init_grass <= 1) 
    init_grass : (float) Initial proportion of grass cover (0 to 1; init_tree + init_grass <= 1)
'''

import numpy as np

def set_probabilities(m): #TODO: J
    ''' set_probilities creates a dictionary of growth, spontaneous ignition, and fire spread probabilities for grass and trees based the input parameters and moisture levels'''
    # Growth rates vary proportionately with moisture level  
        # Set P_ag (probability that grass propagates to neighboring ash site) equal to m 
        # Set P_gt (probability that tree propagates to neighboring grass site) equal to m
    # Ignition probabilities vary inversely with moisture level 
        # Set p_i_t (prob of tree ignition) = p_i_tmax*(1 - m) 
        # Set p_i_g (prob of grass ignition) = p_i_gmax*(1 - m)
    # Fire spread between neighbors depends on moisture and state each cell 
        # Set rate of spread/catch for each cell type equal to (max rate)*(1-m)
        # Multiply (rate of spread) times (rate of catching) for each combination of (grass, tree) to get probabilities of fire spreading between trees and grass 
    # Returns an dict of all the transition probabilities
    p_ag = m
    p_gt = m




    return {"p_ag": p_ag, "p_gt": p_gt}

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


def run_simulation(parameters): #TODO: J
    ''' run_simulation is the wrapper function which initializes the simulation and runs it for a specified number of growth and fire seasons'''
    # Initialize the forest 
    # Set all the probabilities based on moisture levels 
    # For each time steps: 
    # run growth_season() 
    # Run fire_season()
    # Return data 
    pass 


def grow_season(world, params): #TODO: J 
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
    pass


def fire_season(forest, L, params): 
    '''fire_season runs one iteration of the world, in which wildfires initiate, propagate, and die out'''
    # TODO: Make sure there is no bias in the order we iterate through fire cells because
    # we are not descretizing time steps. 

    # -----INITIALIZE----- #
    # get length of non-forest border
    d=int((forest.shape[0] - L)/2)
    # Create two lists: 
    # Fire_location = list of locations of currently on-fire cells in grid
    fire_location=[]
    # Fire_type = 1 or 2 depending on whether each cell is tree or grass fire (aligned with Fire_location list) 
    fire_type=[]

    # -----FIRE SEASON STARTS----- #
    ### IGNITION
    # Randomly select some of the vegetation cells to ignite and add to the lists 
    ignition_probs = np.random.rand(L*L)
    ignition_probs.resize(L,L)
    grass_ignitions = (ignition_probs < params["p_i_g"]) & (forest[d:L+d, d:L+d] == 1)
    tree_ignitions = (ignition_probs < params["p_i_t"]) & (forest[d:L+d, d:L+d] == 2)
    
    # add grass to fire_location list
    # need to add d back in so it lines up with full world with boundaries
    grass_fire_locations = np.where(grass_ignitions) 
    for i, grass_x in enumerate(grass_fire_locations[0]): 
        # Need to add d back in
        fire_location.append((grass_x + d, grass_fire_locations[1][i] + d))
        fire_type.append(1)
    # Add tree to fire_location list
    tree_fire_locations = np.where(tree_ignitions)
    for i, tree_x in enumerate(tree_fire_locations[0]): 
        fire_location.append((tree_x + d, tree_fire_locations[1][i] + d))
        fire_type.append(2)

    ### FIRE SPREADS UNTIL IT IS ALL BURNED OUT
    # While there are still cells on fire:
    while(len(fire_location) > 0):

        # Take the first element in fire_location: 	
        location = fire_location.pop(0)
        type = fire_type.pop(0)
        # set it to ash
        forest[location] = 0

        # Get its neighbors (3x3 moore neighborhood)
        neighbors = forest[location[0]-1:location[0]+2, location[1]-1:location[1]+2]
        # Randomly select which of its neighbors get set on fire (probabilities depend on moisture and focal/neighbor state)
        
        # Get which cells fire spreads to
        # two dimensions, one for spread, one for catch

        # get non-ash neighbors
        foliage_indices = np.nonzero(neighbors)
        foliage_neighbors = neighbors[foliage_indices]
        fire_probs = np.random.rand(2, len(foliage_neighbors))

        if type == 1: # GRASS
            spread = fire_probs[0] < params["r_g_spread"]
        else: # TREE
            spread = fire_probs[0] < params["r_t_spread"]
        
        # Get which cells fire actually ignites, given it spreads
        for i, neighbor in enumerate(foliage_neighbors): 
            if spread[i]: 
                if ( neighbor == 1 # GRASS
                    & (fire_probs[1][i] < params["r_g_catch"]) # grass catch prob
                ) or ( neighbor == 2 # TREE
                      & (fire_probs[1][i] < params["r_t_catch"]) # tree catch prob
                ): 
                    # Get index
                    nx = foliage_indices[0][i]
                    ny = foliage_indices[1][i]
                    new_location = (location[0] + (nx-1), location[1] + (ny-1))
                    # Add to fire list
                    fire_location.append(new_location)
                    # add to type list
                    fire_type.append(neighbor)
                    # change to ash in world

        # (Optionally record data)
    # Return world at next time step (and optional data) 
    return forest
