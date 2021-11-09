"""

This script performs a smoothed inversion on synthetic data to investigate the
ability of the inversion scheme to recover data.

Synthetic input is a 2D sinwave of wavelength specified by user

The set up of the inversion (e.g. topographic, drainage) is the same as used for
the true inversion.

e.g.: python invert_synthetic.py [nx] [ny] [nlam] [eggbox_wavelength]

"""

import netCDF4
import time
import sys
import os
import scipy as sp
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from landlab import RasterModelGrid
from landlab.components import FlowAccumulator, SinkFillerBarnes
from landlab.components.flow_accum.flow_accum_bw import find_drainage_area_and_discharge
from landlab.utils import get_watershed_mask, get_watershed_outlet
from matplotlib.colors import LogNorm

print("Starting at:",datetime.now())
bigstart=time.time()
#### Loading in the input parameters ####

print("Specified inversion parameters:")
print("Number x nodes =", sys.argv[1])
print("Number y nodes =", sys.argv[2])
print("Number of lamda values to try=", sys.argv[3])

#### Set up the model grid with input topography ####

print("Setting up model topographic grid and flowrouting")

zr_nc=netCDF4.Dataset('input_dir/Clyde_Topo_100m_working.nc')
zr_ma = zr_nc['z'][:,:]
mg = RasterModelGrid(zr_ma.shape,xy_spacing=(100,100))

zr = mg.add_zeros('node', 'topographic__elevation')
zr += np.reshape(zr_ma.data.astype(float),zr.shape)
dx = 100

# Set up the boundary conditions on the square grid
mg.set_fixed_value_boundaries_at_grid_edges(True,True,True,True)

flat_shape = zr.shape # a tuple to flatten arrays [number of nodes long]
full_shape = mg.shape # the full shape of the grid [rows, columns]

sfb = SinkFillerBarnes(mg, method='D8',fill_flat=False)
sfb.run_one_step()

'''
plt.figure()
plt.title("Topographic Elevation")
plt.imshow(zr.reshape(full_shape),cmap='cubehelix',origin='lower')
plt.xlabel('Horizontal grid cells')
plt.ylabel('Vertical grid cells')
plt.show()
'''

frr = FlowAccumulator(
    mg,
    'topographic__elevation',
    flow_director = 'FlowDirectorD8')
frr.run_one_step()  # flow routing

a, q_homo_incis = find_drainage_area_and_discharge(mg.at_node['flow__upstream_node_order'], mg.at_node['flow__receiver_node']) # a is number of nodes

area = mg.at_node['drainage_area']
area_threshold = 8
is_drainage = area > (area_threshold*1000000) #km2 to m2

'''
plt.figure()
plt.imshow(is_drainage.reshape(full_shape),origin='lower') 
plt.xlabel('Horizontal grid cells')
plt.ylabel('Vertical grid cells')

plt.title("Channels > 8 km2")
plt.show()
'''

print("Loading in sample localities and geochemical data")

sample_data = np.loadtxt('input_dir/filtered_sample_loc.dat',dtype=str) # [x, y, sample #]
sample_locs = sample_data[:,0:2].astype(float) #get locs as x,y

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

'''
plt.figure()
plt.imshow(is_drainage.reshape(full_shape),origin='lower')
plt.scatter(x=channel_xy[:,0]/mg.dx, y=channel_xy[:,1]/mg.dx,c='grey', s=5)
plt.scatter(x=sample_locs[:,0]/mg.dx, y=sample_locs[:,1]/mg.dx, marker="x",c='r', s=40)
plt.scatter(x=fitted_locs[:,0]/mg.dx, y=fitted_locs[:,1]/mg.dx, marker="+",c='b', s=40)
plt.xlabel('Horizontal grid cells')
plt.ylabel('Vertical grid cells')
plt.title("Manual check of localities (you have to zoom in!)")
plt.show()
'''

loc_indxs = np.transpose(np.flip((fitted_locs/100).astype(int),axis=1))
loc_nodes = np.ravel_multi_index(loc_indxs,dims=full_shape)

obs_data = pd.read_csv('input_dir/converted_chem_data.csv') #read in geochem data
elems =  obs_data.columns[1:].tolist() # List of element strings
obs_data[elems]=obs_data[elems].astype(float) # Cast numeric data to float


print("Setting up inversion grid")

model_width = mg.shape[1] # number of x cells in topographic model grid
model_height = mg.shape[0] # number of y cells in topographic model grid

def expand(block_grid,block_x,block_y):
    """Expands low res array of block heights into
    model grid array that can be fed into topographic
    model. Note that blocks at the upper and eastern
    perimeter are clipped if number of blocks doesn't
    divide number of model cells.

    block_x and block_y are the number of model cells
    in each block in x and y dir respectively"""
    return(block_grid.repeat(block_y, axis=0).repeat(block_x, axis=1)[:model_height,:model_width])

lowest_sample = loc_nodes[sample_data[:,2]== '700012'] # Locality 700012
active_area = get_watershed_mask(mg,lowest_sample) # extract upstream area of most downstream tay sample 

'''
plt.figure()
plt.title("Coverage relative to inversion grid cells")
plt.ylabel("Vertical model grid cells")
plt.xlabel("Horizontal model grid cells")
plt.imshow(active_area.reshape(full_shape),origin='lower')
plt.show()
'''

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

print("Defining synthetic functions")

def eggbox(w,max=1e4,min=1e1):
    """Returns a sin wave chequerboard"""
    out = np.zeros(full_shape)
    for i in np.arange(out.shape[0]):
        for j in np.arange(out.shape[1]):
            out[i,j] += ((((np.sin(i*np.pi/w) * np.sin(j*np.pi/w))+1)/2)*(max-min) +min)
    return(out)

def downstream_synth_data(input_arr):
    """Returns the synthetic data downstream after mixing of input synthetic bedrock array"""
    a, sed_comp = find_drainage_area_and_discharge(mg.at_node['flow__upstream_node_order'], mg.at_node['flow__receiver_node'],runoff = input_arr.reshape(flat_shape))  # composition but homogeneous erosion
    sed_comp_norm = sed_comp/q_homo_incis
    observed = sed_comp_norm[loc_nodes]
    return(observed)

def get_synth_prior(input_arr):
    """Returns the prior that should be used to intitiate the prior based on observations downstream"""
    observed = downstream_synth_data(input_arr)
    clyde_mth_comp =  observed[obs_data['SAMPLE_No'] == 700012]
    synth_prior_wtd_avg = np.nanmean(clyde_mth_comp)
    synth_prior_wtd_avg_log = np.log(synth_prior_wtd_avg)
    return(synth_prior_wtd_avg_log)

print("Defining inverse functions")

def cost_data_synth(bdrck_arr,synthetic_observations):
    """Returns L2norm data misfit for a given bedrock input array (mgkg) relative to synthetic data assuming homogenous mixing"""
    a, sed_comp = find_drainage_area_and_discharge(mg.at_node['flow__upstream_node_order'], mg.at_node['flow__receiver_node'],runoff = bdrck_arr.reshape(flat_shape))  # composition but homogeneous erosion
    sed_comp_norm = sed_comp/q_homo_incis
    l2norm = np.linalg.norm(np.log10(synthetic_observations) - np.log10(sed_comp_norm[loc_nodes]))
    return(l2norm)


# Model roughness
def cost_roughness(blox,active_blox,block_x,block_y):
    """Returns l2 norm of roughness in both directions.
    Assumes Von Neumann BCs dC/dx = dC/dy = 0"""
    copy = np.copy(blox)
    # Set von neumann BCs
    copy[np.invert(active_blox)] = 'nan' # set inactive nodes to 'nan'
    padded = np.pad(copy,pad_width=1,mode='constant',constant_values='nan') # pad with 'nan' too
    x_diffs = np.diff(padded,axis=1)/block_x # dC/dx
    y_diffs = np.diff(padded,axis=0)/block_y # dC/dy
    x_rough = np.sqrt(np.nansum(x_diffs**2)) # sqrt(SUM((dC/dx)^2))
    y_rough = np.sqrt(np.nansum(y_diffs**2)) # sqrt(SUM((dC/dy)^2))
    return(x_rough,y_rough) # return tuple of roughness along both axes

# def cost_model_synth(active_blox,synth_prior):
#     """Returns the l2 norm of the model relative to the prior. input in log space"""
#     l2norm = np.linalg.norm(active_blox - synth_prior)
#     return(l2norm)
#
# cost_model_synth(np.array([1,1,1]),eg_synth_prior)

def synth_smoothed_objective(param_arr,blox,active_blox,block_xstep,block_ystep,synthetic_observations,lamda_):
    """Tests a given parameter array `param_arr` for the given synthetic inversion setup.
    Returns the least squares smoothed cost. Each iteration is ~25 ms"""
    blox[active_blox] = param_arr # Update model grid with new parameters; 1.25 us
    bedrock = expand(np.exp(blox),block_xstep,block_ystep) # Convert log blocks into continuous grid in mg/kg; 3.5 ms
    data_sq = cost_data_synth(bedrock.reshape(flat_shape),synthetic_observations)**2 # Calculate data misfit; 19.4 ms
    rough_x,rough_y = cost_roughness(blox,active_blox,block_xstep,block_ystep) # Calculate roughness; 68 us
    roughness_sq = (lamda_**2)*(rough_x**2 + rough_y**2) # Roughness squared; 0.6 us
    return(data_sq+roughness_sq)

# def synth_damped_objective(param_arr,blox,active_blox,block_xstep,block_ystep,synth_prior,mu_):
#     """Tests a given parameter array `param_arr` for the given synthetic inversion setup.
#     Returns the least squares damped cost. Each iteration is ~25 ms"""
#     blox[active_blox] = param_arr # Update model grid with new parameters; 1.25 us
#     bedrock = expand(np.exp(blox),block_xstep,block_ystep) # Convert log blocks into continuous grid in mg/kg; 3.5 ms
#     data_sq = cost_data_synth(bedrock.reshape(flat_shape),synthetic_observations)**2 # Calculate data misfit; 19.4 ms
#     damping_sq = (cost_model(blox[active_blox],synth_prior)*mu_)**2 # Calculate model size; 15.8 us
#     return(data_sq+damping_sq)

def initiate_synthetic_inversion(nx,ny,prior):
    """Initiates an inversion grid for given number of
    cells and element. """
    blox = np.zeros((ny,nx)) + prior # Set block height to prior in log space
    active_blox = get_active_blocks(nx,ny) # Active cells
    block_x_step = np.ceil(model_width/nx) # Block width
    block_y_step = np.ceil(model_height/ny) # Block height

    return blox,active_blox,block_x_step,block_y_step

print("##### Performing the inversion ######")
eggbox_size=float(sys.argv[4])/dx
# Generate downstream data and priors
source_region = eggbox(eggbox_size)
print(source_region)
synth_obs = downstream_synth_data(source_region)
print(synth_obs)
synth_prior = get_synth_prior(source_region)
'''
plt.imshow(source_region,origin='lower', norm=LogNorm())
plt.show()
'''

# Initiate inversion

nx=int(sys.argv[1])  # <<<<<<<<<<< Change values
ny=int(sys.argv[2])
n_lam = int(sys.argv[3])
lam_exp_list = np.linspace(-2.0,2.0,n_lam)  # given as powers of 10

out_prefix='synth_outputs/smoothed_'+str(nx)+'x'+str(ny)+'_eggbox_'+str(eggbox_size)+'_'

###
#
counter = 0 # <<<<< Uncomment if running as a batch job
array_index= int(os.environ['PBS_ARRAY_INDEX']) # <<<<< Uncomment if running as a batch job

#### Perform inversion ####
for j in np.arange(n_lam):
        counter=counter+1 # <<<<< Uncomment if running as a batch job
        if (counter == array_index): # <<<<< Uncomment if running as a batch job

            blocks,active_blocks,block_width,block_height = initiate_synthetic_inversion(nx,ny,synth_prior)
            parameters = np.copy(blocks[active_blocks]) # The array we will vary.

            #### Perform inversion ####
            lam_exp = lam_exp_list[j]
            lamda = 10**lam_exp

            print("############ Inverting synthetic with lamda = 10^",lam_exp,"############")

            start = time.time()
            print("Starting at: ", datetime.now())
            res_nm = sp.optimize.minimize(fun=synth_smoothed_objective,args=(blocks,active_blocks,block_width,block_height,synth_obs,lamda),x0=parameters,method='Nelder-Mead',
                                      options={'disp':True,'fatol':1e-4,'xatol':1e-4,'maxiter':2e7,'adaptive': True})
            end = time.time()

            #### Finish ####

            print("############ results ############")
            print("Finished at:", datetime.now())
            print("Runtime = ",end-start,"s")
            print(res_nm.success)
            print(res_nm.status)
            print(res_nm.message)
            print(res_nm.nit)

            expanded = expand(np.exp(blocks),block_width,block_height)
            expanded[np.invert(active_area.reshape(full_shape))] = 'nan'

            source_region[np.invert(active_area.reshape(full_shape))] = 'nan'

            blocks[active_blocks] = res_nm.x
            expanded = expand(np.exp(blocks),block_width,block_height)
            expanded[np.invert(active_area.reshape(full_shape))] = 'nan'

            data_misfit = cost_data_synth(expanded.reshape(flat_shape),synth_obs)
            model_roughness = cost_roughness(blocks,active_blocks,block_width,block_height)

            coeff_calib_array = np.array([lam_exp,data_misfit,model_roughness[0],model_roughness[1]])

            # Output interpolated geochemical map
            np.save(out_prefix+str(lam_exp)+'_expanded',expanded)
            # Output node heights
            np.save(out_prefix+str(lam_exp)+'_x',blocks[active_blocks])
            # Output mu calibration array in form: log10(lamda) | datamisfit | roughnessx | roughnessy
            np.savetxt(out_prefix+str(lam_exp)+'_calib.dat',coeff_calib_array.reshape(1, 4))

'''
plt.imshow(expanded,origin='lower',norm=LogNorm())
plt.colorbar()
plt.title("Output Concentrations")
plt.show()
'''