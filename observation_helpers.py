import scipy.ndimage as ndimage
import numpy as np

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

# %%
import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt

ASH = 0 
GRASS = 1
TREE = 2

colors_list = ["black", "yellow", "green"]
cmap = colors.ListedColormap(colors_list)
bounds = [0,1,2,3]
norm = colors.BoundaryNorm(bounds, cmap.N)

forest = np.load('forest1.npy')
plt.imshow(forest,  cmap=cmap, norm=norm)

# Get TREE patch sizes
tree_patches = get_patch_sizes(forest == TREE)

# Get patches of vegetation (TREE or GRASS)
veg_patches = get_patch_sizes(forest == (TREE or GRASS))

# get patches of grass
grass_patches = get_patch_sizes(forest == GRASS)

# %%

