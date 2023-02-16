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

def initialize_forest(L, d, init_grass, init_tree): #TODO: H
    ''' initialize_forest creates the world for the simulation (stored in a 2D array) and populates it with a randomly dispersed initial set of trees and grass'''
    # create (L+2d)x(L+2d) array(forest)
    # Randomly populate the LxL center with trees and grass based on initial conditions
    # Set boundary region (2d length around edge) to be all ash 
    pass 

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


def fire_season(world, params): #TODO: H
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
    pass
