'''
This module gives the code to simulate a forest fire [#TODO: add description]
Written by J Kahn and Helen Miller 
'''
import numpy as np 
from scipy.ndimage import generic_filter


ASH = 0 
GRASS = 1
TREE = 2

# TODO: set seed
# TODO: animate
# TODO: get recurrence interval
# TODO: Save outputs to csv
# TODO: change color scheme
# TODO: find some relevant literature
# TODO: run some parameter sweeps on: 
    # * neighborhood CC (in relationship to tree cc)
    # * growth rate
    # * carrying capacity
    # * ignition probability
    # collect data on mean and sd of biomass (from 200-100 years or something)

def set_probabilities_biomass(p_disp, p_prop, min_seed, r_grow, tree_carrying_capacity, neighborhood_carrying_capacity, max_ignite): 
    '''
    Params: 
        - p_disp: probability of long-range dispersal into any ash cell
        - p_prop: probability that a cell with biomass > min_seed will propogate into neighbor cell
        - min_seed: the minimum biomass that a cell needs before it can propogate
        - r_grow: intrinsic growth rate
        - tree_carrying_capacity: carrying capacity of a cell
        - neighborhood_carrying_capacity: carrying capacity of a neighborhood
        - max_ignite: ignition probability where biomass = 1. Ignition probability = max_ignite / biomass. 
    '''
    return {"p_disp": p_disp, "p_prop": p_prop, "min_seed": min_seed, "r_grow": r_grow, "tree_carrying_capacity": tree_carrying_capacity, "neighborhood_carrying_capacity": neighborhood_carrying_capacity, "max_ignite": max_ignite}


def initialize_params_dict(m, L, t_steps, d, init_grass, init_tree, p_disp, p_prop, min_seed, r_grow, tree_carrying_capacity, neighborhood_carrying_capacity, max_ignite, rand_seed): 
    '''convenience function that puts all the given and calculated parameters into a dictionary, which gets passed 
    around through the simulation''' 
    rng = np.random.default_rng(seed=rand_seed) #initialize a random number generator to create replicable results
    params_dict = {'m': m, 'L' : L, 'd' : d, 'rng': rng}

    probs_dict = set_probabilities_biomass( p_disp=p_disp, 
                                   p_prop=p_prop,
                                   min_seed=min_seed,
                                    r_grow=r_grow, 
                                    tree_carrying_capacity=tree_carrying_capacity,
                                    neighborhood_carrying_capacity=neighborhood_carrying_capacity, 
                                    max_ignite=max_ignite) # Set all the probabilities based on moisture levels 
    params_dict.update(probs_dict) 
    return params_dict

def initialize_forest(L, d, init_grass, init_tree, rng): 
    ''' initialize_forest creates the world for the simulation (stored in a 2D array) and populates it with a randomly dispersed initial set of trees and grass. 
    Ash, grass and trees are represented by ints: 0=ash; 1=grass; 2=tree.'''
    # create (L+2d)x(L+2d) array(forest)
    # Randomly populate forest with trees and grass based on initial conditions
    forest = rng.choice([ASH, GRASS, TREE],
                              size=(L+2*d,L+2*d), 
                              p=[1-(init_tree+init_grass), init_grass, init_tree])

    # Set boundary region (2d length around edge) to be all ash 
    forest[:, 0:d] = ASH
    forest[:, (L+d):] = ASH
    forest[0:d, :] = ASH
    forest[(L+d):, :] = ASH
    return forest    

def count_equals(arr, val): 
    return np.count_nonzero(arr == val)

def grow_season(forest, params_dict):
    '''grow_season runs one iteration of the world in which new vegetation grows'''
    # set up tools + get parameters 
    d = params_dict['d']
    L = params_dict['L']
    p_disp  = params_dict['p_disp']
    p_prop = params_dict['p_prop']
    r_grow = params_dict['r_grow']
    tree_carrying_capacity = params_dict['tree_carrying_capacity']
    neighborhood_carrying_capacity = params_dict['neighborhood_carrying_capacity']
    min_seed = params_dict['min_seed']
    rng = params_dict['rng']

    ### Get all cells where vegetation seeds
    # count sources of maturation age
    
    ### Grow cells with veg!
    # Logistic growth: 
    # individual cells limited by tree carrying capacity. (eg nutrient limits)
    # neighborhood limited by neighborhood carrying capacity (eg light limit)
    # also clip scaling factor at 0,1 so biomass is not lost (even if above carrying capacity)
    # biomass_i(t+1) = r_grow * ( 1 - (biomass_neighborhood(t) / neighborhood_carrying_capacity ) * biomass_i(t)
    # Use 3x3 moore neighborhood for carrying capacity
    neighborhood_biomass = generic_filter(forest, np.sum, size = (3, 3), mode="constant")
    new_grow = r_grow * np.clip(1 - ( (neighborhood_biomass / (neighborhood_carrying_capacity)) + forest / tree_carrying_capacity), 0, 1) * forest

    # probability of grass is probability of long-range dispersal plus probability of 
    # at least 1 neighbor propogating to this square 
    # ie the complement of (NO neighbors propogating there)
    # Initialize biomass at 1. 
    num_propogation_sources = generic_filter(forest > min_seed, np.sum, size = (2*d-1, 2*d-1), mode="constant")
    p_seed = 1 - (1-p_prop) ** num_propogation_sources
    p_seed = np.clip((p_disp + p_seed), 0, 1)
    new_veg = np.where(forest == 0, rng.binomial(1, p_seed), 0)

    forest = forest + new_veg + new_grow

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
    max_ignite = params_dict['max_ignite']
    tree_carrying_capacity = params_dict['tree_carrying_capacity']
    rng = params_dict['rng']
    # Create two lists: 
    # Fire_location = list of locations of currently on-fire cells in grid
    fire_location=[]
    # Fire_type = GRASS or TREE depending on whether each cell is tree or grass fire (aligned with Fire_location list) 
    fire_biomass=[]
    # Keep track of total number of cells burned
    area_burned=0
    # and their locations 
    indices_burned = []

    def ignite(location, biomass): 
        """Helper function which does all the steps needed when a new cell ignites"""
        # Add to fire_location 
        fire_location.append(location)
        # Add to fire_type
        fire_biomass.append(biomass)
        # set location to ash. 
        # technically it's not ash yet, but since we are keeping 
        # track of fire type and location in the other lists we 
        # don't need to in the world. We need to set it to ash 
        # so it doesn't ignite again (ie get doubles in fire_location/fire_type)
        forest[location] = 0

    # -----FIRE SEASON STARTS----- #
    ### IGNITION
    # Randomly select some of the vegetation cells to ignite and add to the lists 
    ignition_probs = rng.random(L*L)
    ignition_probs.resize(L,L)
    # Probability of igniting is max_ignite/biomass
    ignition_probability = max_ignite / forest[d:-d, d:-d]
    # no ignitions at ash sites
    ignition_probability[forest[d:-d, d:-d] == 0] = 0
    ignitions = ignition_probs < ignition_probability
    
    ### end fire season if there were no ignitions 
    if (len(ignitions) == 0): 
        return forest
    
    # add grass ignitions to fire_location list
    fire_locations = np.where(ignitions) 
    for i, x in enumerate(fire_locations[0]): 
        location = (x + d, fire_locations[1][i] + d) # need to add d back in so it lines up with full world with boundaries
        ignite(location=location, biomass = forest[location])

    # -----FIRE SPREADS UNTIL IT IS ALL BURNED OUT----- #
    # While there are still cells on fire:
    while(len(fire_location) > 0):

        # Take the first element in fire_location: 	
        location = fire_location.pop(0)
        biomass = fire_biomass.pop(0)
        area_burned += 1
        indices_burned.append(location)

        # Get its neighbors (3x3 moore neighborhood)
        neighbors = forest[location[0]-1:location[0]+2, location[1]-1:location[1]+2]
        
        # Randomly select which of its neighbors get set on fire (probabilities depend on moisture and focal/neighbor state)
        # get non-ash neighbors
        foliage_indices = np.nonzero(neighbors)
        foliage_neighbors = neighbors[foliage_indices]

        # Probability of ignition = (biomass_source / biomass_target) / tree_carrying_capacity * biomass_source
        
        fire_probs = rng.random(len(foliage_neighbors))

        for i, neighbor in enumerate(foliage_neighbors): 
            p_ignite = (biomass**2 / neighbor) / tree_carrying_capacity
            if p_ignite > fire_probs[i]: 
                # ignite!
                nx = foliage_indices[0][i]
                ny = foliage_indices[1][i]
                # Use relative location to get location in forest
                new_location = (location[0] + (nx-1), location[1] + (ny-1))
                ignite(location=new_location, biomass=neighbor)

                
        # (Optionally record data)
    # Return world at next time step (and optional data) 
    return forest, area_burned, indices_burned


def mask_burned_indices(burn_indices, dim): 
    '''creates and returns a binary mask (1= was burned, 0= not burned)
    Params: 
        - burn_indices = array-like/list of 2-d indices in form [[i1, j1], [i2, j2]...]
        - dim = side length of output mask (should match dim of forest)
    Returns 
        - (nd array) dim x dim mask 
    '''
    mask = np.zeros(dim*dim)
    if len(burn_indices) > 0: # if nothing was burned, return a mask of all 0s
        burn_indices = np.array(burn_indices)
        # convert multi-indices to a list of single dimensional indices 
        burn_area_flat_index = dim*burn_indices[:, 0] + burn_indices[:, 1]
        # apply to a flattened mask
        mask[burn_area_flat_index] = 1
        # reshape mask back to intended dim 
    mask = mask.reshape((dim, dim))
    return mask 


def run_simulation(m, L, t_steps, d, init_grass, init_tree, p_disp, p_prop, min_seed, r_grow, tree_carrying_capacity, neighborhood_carrying_capacity, max_ignite, output_times=[], rand_seed=None): 
    ''' run_simulation is the wrapper function which initializes the simulation and runs it for a specified number of growth and fire seasons
    Params:  
        - m : (float) Moisture level (fundamental parameter dictating probabilities of growth and fire spread) (can vary between 0 and 1, where 0 is no moisture and 1 is total saturation) 
        - L : (int) side length of forest
        - t_steps : (int) Number of “years” (alternations between growth and fire season) to run simulation for
        - d : (int) Radius of neighborhood during growth season (3x3 Moore neighborhood is d=2; exact size of neighborhood given by (2d-1)); also the size of the ash buffer around the edge of the world 
        - p_ig_gmax : (float) Maximum probability of spontaneous grass ignition (when moisture=0) (p_ig_tmax < p_ig_gmax)
        - p_ig_tmax : (float) Maximum probability of spontaneous tree ignition (when moisture=0) (p_ig_tmax < p_ig_gmax)
        - r_spr_tmax : (float) Max rate at which fire spreads from a tree (r_spr_gmax < r_spr_tmax) 
        - r_spr_gmax : (float) Max rate at which fire spreads from grass to neighbors (r_spr_gmax < r_spr_tmax
        - r_cat_tmax : (float) Max rate at which trees catch fire from nearby source (r_cat_tmax < r_cat_gmax)
        - r_cat_gmax : (float) Max rate at which grass catches fire from nearby source (r_cat_tmax < r_cat_gmax)
        - init_tree : (float) Initial proportion of tree cover (0 to 1; init_tree + init_grass <= 1) 
        - init_grass : (float) Initial proportion of grass cover (0 to 1; init_tree + init_grass <= 1)
        - rand_seed : (optional) integer seed for replicable random results 

    Returns:
        - dict with results + data: 
            - output_slices: array, each element contains (alternating) state of forest pre- and post- fire for select years
            - area_burned: array, each element indicates number of vegetation cells burned during each fire season
            - tree_count: array, each element indicates number of tree cells in forest after fire season each year 
            - grass_count: array, each element indicates number of grass cells in forest after fire season each year 
    '''
    params_dict = initialize_params_dict(m=m, 
                                         L=L, 
                                         t_steps=t_steps, 
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
                                    rand_seed=rand_seed)
    forest = initialize_forest(L=L, 
                            d=d, 
                            init_grass=init_grass, 
                            init_tree=init_tree,
                            rng=params_dict['rng'])

    save_counter = 0 # save state data at select times 
    output_slices = np.zeros((len(2*output_times), (L+2*d), (L+2*d)))

    # initialize some key things to save at every time point
    area_burned_mask = np.zeros((t_steps, (L+2*d), (L+2*d)))
    area_burned_output = np.zeros(t_steps)
    total_biomass = np.zeros(t_steps)
    for t in range(t_steps): 
        forest = grow_season(forest, params_dict) #TODO: return data?

        # save select outputs pre-fire 
        if t in output_times: 
            output_slices[save_counter] = forest
            save_counter += 1

        forest, area_burned, indices_burned = fire_season(forest, params_dict)

        # todo: make indices_burned into an np array mask that can be applied over the state
        # Save outputs which get recorded every year
        area_burned_output[t] = area_burned
        total_biomass[t] = np.sum(forest)

        # Save select outputs post-fire 
        if t in output_times: 
            output_slices[save_counter] = forest
            area_burned_mask[save_counter] = mask_burned_indices(indices_burned, L+2*d)
            save_counter += 1
            
    return {"output_times": output_times, 
            "output_slices": output_slices, 
            "burned_area_masks" : area_burned_mask,
            "area_burned": area_burned_output, 
            "biomass": total_biomass,
            "params_dict": params_dict}
