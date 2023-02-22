'''
This module gives the code to simulate a forest fire [add description]
Written by J Kahn and Helen Miller 
'''
import numpy as np 

ASH = 0 
GRASS = 1
TREE = 2

def set_probabilities(m, p_ig_gmax, p_ig_tmax, r_spr_tmax, r_spr_gmax, r_cat_tmax, r_cat_gmax): 
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

def initialize_forest(L, d, init_grass, init_tree): 
    ''' initialize_forest creates the world for the simulation (stored in a 2D array) and populates it with a randomly dispersed initial set of trees and grass. 
    Ash, grass and trees are represented by ints: 0=ash; 1=grass; 2=tree.'''
    # create (L+2d)x(L+2d) array(forest)
    # Randomly populate forest with trees and grass based on initial conditions
    forest = np.random.choice([ASH, GRASS, TREE],
                              size=(L+2*d,L+2*d), 
                              p=[1-(init_tree+init_grass), init_grass, init_tree])

    # Set boundary region (2d length around edge) to be all ash 
    forest[:, 0:d] = ASH
    forest[:, (L+d):] = ASH
    forest[0:d, :] = ASH
    forest[(L+d):, :] = ASH
    return forest

def get_neighbors(forest, location, r):
    '''Retrieves the values of the neighbors of the cell at the given location 
    in the forest

    Params: 
        - forest: 2d array
        - location: (row, col) tuple
        - r : radius of neighborhood, r>=1 (for just the focal cell, r=1; for 3x3 Moore neighborhood, r=2)
    Returns: 
        - 2d array with dimensions (2r-1, 2r-1)
    '''
    low_bound = r-1 
    up_bound = r
    return forest[location[0]-low_bound:location[0]+up_bound, location[1]-low_bound:location[1]+up_bound]

def try_propagate(neighbors, p_prop, current_val, new_val):
    """
    attempts to propagate vegetation from surrounding cells into the focal cell. returns 
    the new state of the focal cell 
    Args:
        neighbors (nd_array): array giving all the neighbors of the focal cell 
        p_prop (float): probability of new vegetation growing in the focal cell from 1 neighbor 
        current_val (int): the current state of the focal cell 
        new_val (int): the 'spreading' vegetation

    Returns:
        int: the value of the cell after the propagation attempt
    
    probability of new veg growing in this square is equal to the probability 
    of at least 1 neighbor propagating to this square 
    ie the complement of (NO neighbors propagating here)
    """
    rng = np.random.default_rng()
    num_surrounding_sources = np.count_nonzero(neighbors == new_val)
    p_grow = 1 - (1-p_prop)**(num_surrounding_sources)
    if rng.binomial(n=1, p=p_grow): # random trial to see if new plant grows
        return new_val 
    else: 
        return current_val


def grow_season(forest, params_dict):
    '''grow_season runs one iteration of the world in which new vegetation grows'''
    # set up tools + get parameters 
    forest_iter = np.nditer(forest, flags=['multi_index']) # BUG: does not respect border boundary conditions 
    d = params_dict['d']
    L = params_dict['L']
    p_gro_ag = params_dict['p_gro_ag'] 
    p_gro_gt = params_dict['p_gro_gt']

    # Iterate across the forest grid 
    for forest_cell in forest_iter: 
        current_index = forest_iter.multi_index
        neighbors = get_neighbors(forest, current_index, d) #get states of neighbors in (2d-1)x(2d-1) neighborhood
       
        # If ash: roll to see if grows grass or stays ash
        if forest_cell == ASH: 
          forest[current_index] = try_propagate(neighbors, p_gro_ag, ASH, GRASS)
        # If grass: roll to see if grows tree or stays grass 
        if forest_cell == GRASS: 
            forest[current_index] = try_propagate(neighbors, p_gro_gt, GRASS, TREE) 

    # make the buffer stay ash
    forest[0:d, ] = ASH # top strip
    forest[:, 0:d] = ASH # left strip
    forest[-d:, ] = ASH # bottom strip
    forest[:, -d: ] = ASH # right strip 
    # Return world at next time step #TODO: collect data 
    return forest

def fire_season(forest, params_dict): 
    '''fire_season runs one iteration of the world, in which wildfires initiate, propagate, and die out'''
    # TODO: add some checks/tests so function throws useful errors if something goes wrong
    
    # -----INITIALIZE----- #
    L=params_dict["L"]
    d=params_dict["d"]
    # Create two lists: 
    # Fire_location = list of locations of currently on-fire cells in grid
    fire_location=[]
    # Fire_type = 1 or 2 depending on whether each cell is tree or grass fire (aligned with Fire_location list) 
    fire_type=[]
    # Keep track of total number of cells burned
    area_burned=0

    def ignite(location, veg_type): 
        """Helper function which does all the steps needed when a new cell ignites"""
        # Add to fire_location 
        fire_location.append(location)
        # Add to fire_type
        fire_type.append(veg_type)
        # set location to ash. 
        # technically it's not ash yet, but since we are keeping 
        # track of fire type and location in the other lists we 
        # don't need to in the world. We need to set it to ash 
        # so it doesn't ignite again (ie get doubles in fire_location/fire_type)
        forest[location] = ASH

    # -----FIRE SEASON STARTS----- #
    ### IGNITION
    # Randomly select some of the vegetation cells to ignite and add to the lists 
    ignition_probs = np.random.rand(L*L)
    ignition_probs.resize(L,L)
    grass_ignitions = (ignition_probs < params_dict["p_ig_g"]) & (forest[d:L+d, d:L+d] == GRASS)
    tree_ignitions = (ignition_probs < params_dict["p_ig_t"]) & (forest[d:L+d, d:L+d] == TREE)
    
    if (len(grass_ignitions)+len(tree_ignitions) == 0): 
        return forest
    
    # add grass to fire_location list
    # need to add d back in so it lines up with full world with boundaries
    grass_fire_locations = np.where(grass_ignitions) 
    for i, grass_x in enumerate(grass_fire_locations[0]): 
        # Need to add d back in
        location = (grass_x + d, grass_fire_locations[1][i] + d)
        ignite(location=location, veg_type=1)
    # Add tree to fire_location list
    tree_fire_locations = np.where(tree_ignitions)
    for i, tree_x in enumerate(tree_fire_locations[0]): 
        location = (tree_x + d, tree_fire_locations[1][i] + d)
        ignite(location=location, veg_type=2)

    ### FIRE SPREADS UNTIL IT IS ALL BURNED OUT
    # While there are still cells on fire:
    while(len(fire_location) > 0):

        # Take the first element in fire_location: 	
        location = fire_location.pop(0)
        veg_type = fire_type.pop(0)
        area_burned += 1

        # Get its neighbors (3x3 moore neighborhood)
        neighbors = forest[location[0]-1:location[0]+2, location[1]-1:location[1]+2]
        
        # Randomly select which of its neighbors get set on fire (probabilities depend on moisture and focal/neighbor state)
        # get non-ash neighbors
        foliage_indices = np.nonzero(neighbors)
        foliage_neighbors = neighbors[foliage_indices]

        ### Using spread and catch rates
        if False: 
            # Get which cells fire spreads to
            # two dimensions, one for spread, one for catch
            fire_probs = np.random.rand(2, len(foliage_neighbors))

            if veg_type == GRASS: 
                spread = fire_probs[0] < params_dict["r_g_spread"]
            else: # TREE
                spread = fire_probs[0] < params_dict["r_t_spread"]
            
            # Get which cells fire actually ignites, given it spreads
            for i, neighbor in enumerate(foliage_neighbors): 
                if spread[i]: 
                    if ( neighbor == GRASS
                        & (fire_probs[1][i] < params_dict["r_g_catch"]) # grass catch prob
                    ) or ( neighbor == TREE
                        & (fire_probs[1][i] < params_dict["r_t_catch"]) # tree catch prob
                    ): 
                        # Get index
                        nx = foliage_indices[0][i]
                        ny = foliage_indices[1][i]
                        new_location = (location[0] + (nx-1), location[1] + (ny-1))
                        ignite(location=new_location, type=neighbor)
        
        # with p_spr_*: 
        if True: 
            fire_probs = np.random.rand(len(foliage_neighbors))

            for i, neighbor in enumerate(foliage_neighbors): 
                if (veg_type == GRASS
                        and (neighbor == GRASS
                            and (fire_probs[i] < params_dict["p_spr_gg"])) 
                        or  (neighbor == TREE 
                            and (fire_probs[i] < params_dict["p_spr_gt"]))
                    ) or (veg_type == TREE
                        and (neighbor == GRASS
                            and (fire_probs[i] < params_dict["p_spr_tg"])) 
                        or  (neighbor == TREE
                            and (fire_probs[i] < params_dict["p_spr_tt"]))
                    ): 
                    
                    # ignite!
                    nx = foliage_indices[0][i]
                    ny = foliage_indices[1][i]
                    # Use relative location to get location in forest
                    new_location = (location[0] + (nx-1), location[1] + (ny-1))
                    ignite(location=new_location, veg_type=neighbor)
                      


        # (Optionally record data)
    # Return world at next time step (and optional data) 
    return forest, area_burned

def initialize_params_dict(m, L, t_steps, d, init_grass, init_tree, p_ig_gmax, p_ig_tmax, r_spr_tmax, r_spr_gmax, r_cat_tmax, r_cat_gmax): 
    '''convenience function that puts all the given and calculated parameters into a dictionary, which gets passed 
    around through the simulation''' 
    params_dict = {'m': m, 'L' : L, 'd' : d}
    probs_dict = set_probabilities(m, p_ig_gmax, p_ig_tmax, r_spr_tmax, r_spr_gmax, r_cat_tmax, r_cat_gmax) # Set all the probabilities based on moisture levels 
    params_dict.update(probs_dict) 
    return params_dict

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
    params_dict = initialize_params_dict(m, L, t_steps, d, init_grass, init_tree, p_ig_gmax, p_ig_tmax, r_spr_tmax, r_spr_gmax, r_cat_tmax, r_cat_gmax)
    
    save_counter = 0 # save state data at select times 
    output_slices = np.zeros((len(output_times), (L+2*d), (L+2*d)))

    # initialize some key things to save at every time point
    area_burned_output = np.zeros(t_steps)
    tree_count = np.zeros(t_steps)
    grass_count = np.zeros(t_steps)
    for t in range(t_steps): 
        forest = grow_season(forest, params_dict) #TODO: return data 
        forest, area_burned = fire_season(forest, params_dict)

        # Save outputs which get recorded every year
        area_burned_output[t] = area_burned
        tree_count[t] = np.sum(forest == TREE)
        grass_count[t] = np.sum(forest == GRASS)

        # Save select outputs
        if t in output_times: 
            output_slices[save_counter] = forest
            save_counter += 1
            
    return {"output_slices": output_slices, 
            "area_burned": area_burned_output, 
            "tree_count": tree_count, 
            "grass_count": grass_count}

