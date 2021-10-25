import netCDF4
import sys
import time
import scipy as sp
import landlab
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from landlab import RasterModelGrid
from landlab.components import FlowAccumulator, FastscapeEroder, SinkFillerBarnes, FlowDirectorD8, ChannelProfiler, TrickleDownProfiler
from landlab.components.flow_accum.flow_accum_bw import find_drainage_area_and_discharge
from landlab import imshow_grid
from landlab.utils import get_watershed_mask,get_watershed_outlet,get_watershed_nodes
import argparse
import os

nb_output = sys.stdout # Location to write to console not the notebook
console_output = open('/dev/stdout', 'w') # Location to write to console

################################################## command line parser #########################################
parser = argparse.ArgumentParser()
parser.add_argument('Element', type=str, help= 'Element to invert for, must be given in the same way as it is found in the observations dataframe')
parser.add_argument('minimum_lambda', type=float, help='Minimum minimum lambda (smoothing factor) to try')
parser.add_argument('maximum_lambda', type=float, help='Maximum exponent of lambda (smoothing factor) to try')
parser.add_argument('number_lambda', type=int, help= 'Number of lambda values to try out')

args = parser.parse_args()
elem = args.Element
min_lambda = args.minimum_lambda
max_lambda = args.maximum_lambda
num_lambda = args.number_lambda

################################################## model function definitions ##################################
def expand(block_grid,block_x,block_y):
    """Expands low res array of block heights into 
    model grid array that can be fed into topographic
    model. Note that blocks at the upper and eastern 
    perimeter are clipped if number of blocks doesn't 
    divide number of model cells. 
    
    block_x and block_y are the number of model cells 
    in each block in x and y dir respectively"""
    return(block_grid.repeat(block_y, axis=0).repeat(block_x, axis=1)[:model_height,:model_width])

def get_active_blocks(nx,ny):
    """For a given number of blocks in the x 
    and y direction (nx & ny), returns a (ny,nx) 
    bool array saying if cell overlaps with active 
    area or not. """
    
    block_width = np.ceil(model_width/nx) 
    block_height = np.ceil(model_height/ny)

    model_grid_block_indices = np.zeros((model_height,model_width,2))    
    for i in np.arange(model_height):
        for j in np.arange(model_width):
            model_grid_block_indices[i,j,0] = i//block_height
            model_grid_block_indices[i,j,1] = j//block_width
    model_grid_block_indices = model_grid_block_indices.astype(int)        
    # 3D array that contains index of block that model cell corresponds to 
    # ([:,:,0] = y coordinate; [:,:,1] = x coordinate)
    
    out = np.zeros((ny,nx)).astype(bool)
    for i in np.arange(ny):
        for j in np.arange(nx):
            # Boolean array of model cells that correspond to block indeix (i,j)
            cells_in_block = np.logical_and(model_grid_block_indices[:,:,0] == i, model_grid_block_indices[:,:,1] == j)
            # Returns if block overlap with active area:
            out[i,j] = np.any(np.logical_and(cells_in_block,active_area.reshape(full_shape)))
    return(out)
def loc_to_area(loc):
    """Returns the catchment mask for a given locality"""
    return(loc_areas[np.where(unique_locs == loc)])

def b_is_nested_in_a(a,b):
    """Is catchment 'b' a nested subcatchment of catchment 'a'?"""
    return(not(np.any(np.invert(np.logical_or(a,np.invert(b))))))

def which_locs_are_contained_upstream_of(loc_a):
    """Which localities define subcatchments of area defined by loc 'a'"""
    loc_area = loc_to_area(loc_a)
    out_locs = []
    for loc_num in unique_locs:
        if(not(loc_num==loc_a)):
            upst_area = loc_to_area(loc_num)
            if(b_is_nested_in_a(loc_area,upst_area)):
                out_locs = out_locs + [loc_num]
    return(out_locs)

def find_unique_seg(loc):
    locs_upstream_of_ = which_locs_are_contained_upstream_of(loc)
    downstream_area = loc_to_area(loc)
    out = np.zeros(active_area.size).astype(bool)
    for upstream_loc in locs_upstream_of_:
        out = np.logical_or(out,loc_to_area(upstream_loc))
    unique_seg = np.logical_and(downstream_area,np.invert(out))
    return(unique_seg.reshape(active_area.shape))

def data_misfit(bdrck_arr,elem):
    """Returns L2norm data misfit for a given bedrock input array (mgkg), calculating predictions for given element assuming homogenous incision"""
    a, sed_comp = find_drainage_area_and_discharge(mg.at_node['flow__upstream_node_order'], mg.at_node['flow__receiver_node'],runoff = bdrck_arr)  # composition but homogeneous erosion
    sed_comp_norm = sed_comp/q_homo_incis
    l2norm = np.linalg.norm(np.log10(obs_data[elem]) - np.log10(sed_comp_norm[loc_nodes]))
    return(l2norm)   

def cost_roughness(blox,active_blox,block_x,block_y):
    """Returns l2 norm of roughness in both directions. 
    Assumes Von Neumann BCs dC/dx = dC/dy = 0"""
    copy = np.copy(blox)
    # Set von neumann BCs
    copy[np.invert(active_blox)] = 'nan' # set inactive nodes to 'nan'
    padded = np.pad(copy,pad_width=1,mode='constant',constant_values='nan') # pad with 'nan' too
    x_diffs = np.diff(padded,axis=1)/block_x # dC/dx
    y_diffs = np.diff(padded,axis=0)/block_y # dC/dy
    x_rough = np.sqrt(np.nansum(x_diffs**2)) # sqrt(SUM((dC/dx)^2)), NB 'nans' are treated as zeros using nansum.
    y_rough = np.sqrt(np.nansum(y_diffs**2)) # sqrt(SUM((dC/dy)^2))
    return(x_rough,y_rough) # return tuple of roughness along both axes

def smoothed_objective(param_arr,blox,active_blox,block_xstep,block_ystep,elem,lamda_):
    """Tests a given parameter array `param_arr` for the given inversion setup. 
    Returns the least squares damped cost. Each iteration is ~25 ms"""
    blox[active_blox] = param_arr # Update model grid with new parameters; 1.25 us
    bedrock = expand(np.exp(blox),block_xstep,block_ystep) # Convert log blocks into continuous grid in mg/kg; 3.5 ms
    data_sq = data_misfit(bedrock.reshape(flat_shape),elem)**2 # Calculate data misfit; 19.4 ms
    rough_x,rough_y = cost_roughness(blox,active_blox,block_xstep,block_ystep) # Calculate roughness; 68 us
    roughness_sq = (lamda_**2)*(rough_x**2 + rough_y**2) # Roughness squared; 0.6 us
    return(data_sq+roughness_sq)

def initiate_blocky_inversion_smart(nx,ny,elem):
    """Initiates an inversion grid for given number of 
    cells and element. """
    
    # Define full-res starting solution
    full_init = np.zeros(active_area.shape) + prior_wtd_avg_log[elem]  
    for loc_num in unique_locs:
        values = np.asarray(obs_elems[elem][obs_data['SAMPLE_No'] == float(loc_num)])
        if(values.size==2): # Catch duplicate sample exception
            full_init[find_unique_seg(loc_num)] = np.mean(np.log(values))
        else:
            full_init[find_unique_seg(loc_num)] = np.log(values)

    # Define inversion nodes
    blox = np.zeros((ny,nx))
    active_blox = get_active_blocks(nx,ny) # Active cells
    
    block_x_step = np.ceil(model_width/nx) # Block width
    block_y_step = np.ceil(model_height/ny) # Block height
    
    # Downsample initial guess
    model_grid_block_indices = np.zeros((model_height,model_width,2))
    for i in np.arange(model_height):
        for j in np.arange(model_width):
            model_grid_block_indices[i,j,0] = i//block_y_step
            model_grid_block_indices[i,j,1] = j//block_x_step
    model_grid_block_indices = model_grid_block_indices.astype(int)
    # 3D array that contains index of block that model cell corresponds to
    # ([:,:,0] = y coordinate; [:,:,1] = x coozzzrdinate)
     
    for i in np.arange(ny):
        for j in np.arange(nx):
            # Boolean array of model cells that correspond to block indeix (i,j)
            cells_in_block = np.logical_and(model_grid_block_indices[:,:,0] == i, model_grid_block_indices[:,:,1] == j)
            # Returns if block overlap with active area:
            blox[i,j] = np.mean(full_init.reshape(full_shape)[cells_in_block])
    blox[np.invert(active_blox)] =  prior_wtd_avg_log[elem]     
    return blox,active_blox,block_x_step,block_y_step


########################################### setting up the topographic data ######################################
zr_nc=netCDF4.Dataset('input_dir/Clyde_Topo_100m_working.nc')
zr_ma = zr_nc['z'][:,:]
plt.imshow(zr_ma, origin='lower')
mg = RasterModelGrid(zr_ma.shape,xy_spacing=(100,100))

zr = mg.add_zeros('node', 'topographic__elevation')
zr += np.reshape(zr_ma.data.astype(float),zr.shape)
dx = 100

# Set up the boundary conditions on the square grid
mg.set_fixed_value_boundaries_at_grid_edges(True,True,True,True)

flat_shape = zr.shape # a tuple to flatten arrays [number of nodes long]
full_shape = mg.shape # the full shape of the grid [rows, columns]

#sink filling algorithm: only needs to be run once, otherwise load saved elevationfrom cell below
sfb = SinkFillerBarnes(mg, method='D8',fill_flat=False)
sfb.run_one_step()

#extract data from mg to numpy array (at node)
elevation = mg.field_values('node', 'topographic__elevation')
elevation = elevation.reshape(mg.shape)

# instantiate the flow routing:
frr = FlowAccumulator(
    mg,
    'topographic__elevation',
    flow_director = 'FlowDirectorD8')
frr.run_one_step()  # flow routing

zrold = zr[mg.nodes] # pre-incision elevations

a, q_homo_incis = find_drainage_area_and_discharge(mg.at_node['flow__upstream_node_order'], mg.at_node['flow__receiver_node']) # a is number of nodes
# q is cumulative flux for homogeneous incision

area = mg.at_node['drainage_area']
np.amax(area)
area_threshold = 8 #float(sys.argv[1]) #25 km2
is_drainage = area > (area_threshold*1000000) #km2 to m2
mg.add_field('node','channels',is_drainage,clobber=True)

###################################### loading in chemistry data and snapping to grid ################################
sample_data = np.loadtxt('input_dir/filtered_sample_loc.dat',dtype=str) # [x, y, sample #]
sample_locs = sample_data[:,0:2].astype(np.float)
channel_xy = np.flip(np.transpose(np.where(is_drainage.reshape(mg.shape))),axis=1)*100 # xy coordinates of channels
nudge = np.zeros(sample_locs.shape) # initiate nudge array

#nudging locations to snap to correct channel
nudge[60] = [0,-400]    #nudging loc 700000 to S
nudge[17] = [0,-200]    #nudging loc 632137 to S
nudge[34] = [-700,0]    #nudging loc 632164 to W
nudge[38] = [0,-400]    #nudging loc  632170 to S
nudge[39] = [-100,0]    #nudging loc 632171 to W
nudge[56] = [0,100]     #nudging loc 632197 to N
nudge[16] = [-300,-100] #nudging loc 632136 to SW
nudge[4 ] = [-300,-100] #nudging loc 632109 to SW
nudge[50] = [0,-100]    #nudging loc 632189 to S
nudge[3 ] = [-200,-100] #nudging loc 632108 to SW
nudge[70] = [0,100]     #nudging loc 700012 to N
nudge[66] = [0, -100]

nudged_locs = sample_locs + nudge # Apply the nudges
# Fit the data to the nearest channel node
fitted_locs = np.zeros(sample_locs.shape) # Initialise snapped locality array
for i in np.arange(nudged_locs.shape[0]):
    sample = nudged_locs[i,:]
    diffs = channel_xy-sample
    distances = np.linalg.norm(diffs,axis=1)
    shortest = np.amin(distances)
    fitted_locs[i,:] = channel_xy[np.where(distances==shortest)]

loc_indxs = np.transpose(np.flip((fitted_locs/100).astype(int),axis=1))
print(np.sum(fitted_locs%100))
loc_nodes = np.ravel_multi_index(loc_indxs,dims=full_shape)

# Following statement should be true if all localities are correctly on model grid
print((is_drainage[loc_nodes]).all())

# loc_nodes is used to quick access data in the model grid at sample localities. 
# For example, we now access the *locality* number of every sample via the grid

locality_num_grid = np.zeros(flat_shape) - 20
locality_num_grid[loc_nodes] = sample_data[:,2].astype(float).astype(int)

mg.add_field('node','loc_nums',locality_num_grid,clobber=True)

obs_data = pd.read_csv('input_dir/converted_chem_data.csv') #read in geochem data
elems =  obs_data.columns[1:].tolist() # List of element strings
obs_data[elems]=obs_data[elems].astype(float) # Cast numeric data to float

#loading in compositional data:
obs_elems = obs_data[elems]
clyde_mth_comp = np.asarray(obs_elems[obs_data['SAMPLE_No'] == 700012])


prior_wtd_avg = pd.DataFrame(clyde_mth_comp)
prior_wtd_avg.columns = elems

prior_wtd_avg = np.mean(prior_wtd_avg,axis=0)
prior_wtd_avg_log = np.log(prior_wtd_avg)

#takes long
unique_locs = np.unique(sample_data[:,2])

loc_areas = []
for loc_num in unique_locs:
    print(loc_num)
    sample_node_num = loc_nodes[sample_data[:,2]==loc_num]
    if(sample_node_num.size==2): # Catch duplicate sample exception
        sample_node_num = sample_node_num[0]
    upst_area = get_watershed_mask(mg,sample_node_num)
    loc_areas = loc_areas + [upst_area]
loc_areas = np.array(loc_areas) # The full (not unique) upstream area for each sample site.

model_width = mg.shape[1] # number of x cells in topographic model grid
model_height = mg.shape[0] # number of y cells in topographic model grid
lowest_sample = loc_nodes[sample_data[:,2]== '700012'] # Locality 700012
active_area = get_watershed_mask(mg,lowest_sample) # extract upstream area of most downstream tay sample 
############################# running the model looping through to find best smoothing ######################

nx=10 #  # <<<<<<<<<<<<<<<<   Change number of x blocks in inversion grid
ny=10 #  # <<<<<<<<<<<<<<<<   Change number of y blocks in inversion grid

# crate result array: 2D array [[lambda,rms]]

lambda_exp_array = np.linspace(min_lambda,max_lambda,num_lambda) # Change smoothing factors to try
blocks,active_blocks,block_width,block_height = initiate_blocky_inversion_smart(nx,ny,elem) #initiate inversion blocks
counter = 0 #initiate counter for batch job
array_index = int(os.environ['PBS_ARRAY_INDEX'])

for l in lambda_exp_array:

    counter += 1
    if(counter == array_index): #run single inversion on each node 
        ### Initiating inversion ####
        parameters = np.copy(blocks[active_blocks]) # The array we will vary. 
        worked_lambda = 10**l #raise l to power of 10 to get real lambda value
        worked_exp_lambda = l #store l used in inversion on a given node

        #### Perform inversion ####
        start = datetime.now()
        res_nm = sp.optimize.minimize(fun=smoothed_objective,args=(blocks,active_blocks,block_width,block_height,elem,worked_lambda),x0=parameters,method='Powell',
                                    options={'disp':True,'xtol':1e-3,'ftol':1e-3}) #raising l to power of 10
        end = datetime.now()
        delta=end-start
        print("Total time taken: %.2f" %delta.total_seconds())

        #### Finish ####
        expanded = expand(np.exp(blocks),block_width,block_height)
        expanded[np.invert(active_area.reshape(full_shape))] = 'nan'

        #save result as .asc file with name elem_lambda_inverse_result.asc
        output_path = 'inverse_results/' + elem + '_' + worked_exp_lambda +'_inverse_output.asc'
        np.savetxt(output_path,expanded)

x_rough, y_rough  = cost_roughness(blocks, active_blocks,block_width,block_height) #calculating roughness of output
final_misfit = data_misfit(expanded, elem) #calculating misfit of output
tradeoff_array = np.array(worked_lambda, x_rough, y_rough,final_misfit)
tradeoff_path = 'inverse_results/' + elem + '_' + worked_exp_lambda + '_roughness_misfit.txt'
np.savetxt(tradeoff_path,tradeoff_array, fmt=[ '%d', '%.5f', '%.5f', '%.5f']) #saving roughness and misfit to txt file for later concatination and use