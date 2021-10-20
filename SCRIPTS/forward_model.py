import ipywidgets as widgets
import netCDF4
import sys
import time
from numpy.lib.function_base import diff
import scipy as sp
import landlab
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from landlab import RasterModelGrid
from landlab.components import FlowAccumulator, FastscapeEroder, SinkFillerBarnes, FlowDirectorD8, ChannelProfiler, TrickleDownProfiler
from landlab.components.flow_accum.flow_accum_bw import find_drainage_area_and_discharge
from landlab import imshow_grid
from landlab.utils import get_watershed_mask,get_watershed_outlet,get_watershed_nodes
nb_output = sys.stdout # Location to write to console not the notebook
console_output = open('/dev/stdout', 'w') # Location to write to console

###########################################INPUTS###############################################
element='Mg' #<<<<<<<<<<<<<<<<<<<<<<<< change to correct element, REMEMBER TO INTERPOLATE FIRST
interpolate_input = 'msci-clyde-unmixing/DATA/INTERPOLATED_GBASE/gbase_log_' + element + '.nc' #path to interpolated G-BASE data
result_output = 'msci-clyde-unmixing/DATA/FORWARDMODEL_RESULTS/gbase_log_sed' + element + '.asc' #path to full saved output
misfit_output = 'msci-clyde-unmixing/DATA/FORWARDMODEL_RESULTS/' + element + '_obs_v_pred.txt' #path to output at observed localities


#loading in filled topographic data:
zr_nc=netCDF4.Dataset('msci-clyde-unmixing/DATA/Clyde_Topo_100m_working.nc')
zr_ma = zr_nc['z'][:,:]
topography = np.load('msci-clyde-unmixing/SCRIPTS/filled_topography.npy')

mg = RasterModelGrid(zr_ma.shape,xy_spacing=(100,100)) #add array as topography field
zr = mg.add_field('node', 'topographic__elevation', topography)
dx = 100
flat_shape = zr.shape # a tuple to flatten arrays [number of nodes long]
full_shape = mg.shape # the full shape of the grid [rows, columns]

# Set up the boundary conditions on the square grid
mg.set_fixed_value_boundaries_at_grid_edges(True,True,True,True)

#run flow routing algorithm:
frr = FlowAccumulator(
    mg,
    'topographic__elevation',
    flow_director = 'FlowDirectorD8')
frr.run_one_step()  # flow routing

#Find channels with minimum drainage area of 8km^2 and calculate incision
area = mg.at_node['drainage_area']
np.amax(area)
area_threshold = 8 #float(sys.argv[1]) #25 km2
is_drainage = area > (area_threshold*1000000) #km2 to m2
mg.add_field('node','channels',is_drainage,clobber=True)

#load in interpolated G-BASE data:
comp_nc=netCDF4.Dataset(interpolate_input)
comp_ma = comp_nc['z'][:,:]
comp_log = mg.add_zeros('node','log_bdrck')
comp_log += np.reshape((comp_ma.data.astype(float)),comp_log.shape)
comp = mg.add_zeros('node','bdrck')
comp += np.reshape(((10**comp_ma).data.astype(float)),comp.shape) # Convert to raw values from log10

#Find total sediment flux first, assuming homogenous incision:
a, q = find_drainage_area_and_discharge(mg.at_node['flow__upstream_node_order'], mg.at_node['flow__receiver_node']) # a is number of nodes


#Run forward model using composition and homogenous erosion
a, sed_comp = find_drainage_area_and_discharge(mg.at_node['flow__upstream_node_order'], mg.at_node['flow__receiver_node'],runoff = comp)
sed_norm = sed_comp/q #normalise composition by total sediment flux
sed_norm[q==0] = comp[q==0] #setting composition to bedrock composition where sed flux is 0


#visualise by turning back to log10 and running through channel system:
sed_comp_norm_channel = np.log10(sed_norm) * is_drainage

print(np.any(np.isnan(sed_comp_norm_channel)))
#Add model results to raster model grid:
mg.add_field('node','homo_incis_sed',sed_norm,noclobber=False)
mg.add_field('node','homo_incis_log_sed',np.log10(sed_norm),noclobber=False)
mg.add_field('node','homo_incis_log_sed_channel',sed_comp_norm_channel,noclobber=False)

#saving the result:
#mg.save(result_output, names=['homo_incis_log_sed_channel'])

####################################calculating data misfit between observations and predictions:######################
sample_data = np.loadtxt('msci-clyde-unmixing/DATA/filtered_sample_loc.dat',dtype=str) # [x, y, sample #]
sample_locs = sample_data[:,0:2].astype(np.float)
channel_xy = np.flip(np.transpose(np.where(is_drainage.reshape(mg.shape))),axis=1)*100 # xy coordinates of channels
nudge = np.zeros(sample_locs.shape) # initiate nudge array

#nudging locations:
#indices of locs to nudge: 3,4,16,17,34,38,39,50,50,56,60, 70
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

obs_data = pd.read_csv('msci-clyde-unmixing/DATA/converted_chem_data.csv') #read in geochem data
elems =  obs_data.columns[1:].tolist() # List of element strings
obs_data[elems]=obs_data[elems].astype(float) # Cast numeric data to float

#creating arrays of observations and predictions
obs_log = np.log10(obs_data[element])
pred_log = sed_comp_norm_channel[loc_nodes]
diff_array = obs_log-pred_log
print(np.sqrt(np.sum(diff_array)**2/len(obs_log))) #printint rms
out_array = np.array([sample_data[:,2].astype(int), sample_data[:,0].astype(float), sample_data[:,1].astype(float), pred_log, obs_log, diff_array]).T #create output array that has the form sample#, x, y, obs, pred, obs-pred
np.savetxt(misfit_output, out_array, fmt = ['%d', '%.18f', '%.18f','%.5f', '%.5f', '%.5f'])

#create river profile:
node_eg = mg.grid_coords_to_node_id(377,769) #random node inside catchment
outlet = get_watershed_outlet(mg,node_eg) #finding sink node of Clyde

profiler = ChannelProfiler(mg,main_channel_only=True,outlet_nodes=[outlet]) # Initiates the long profiler 
profiler.run_one_step() # Extract longest profile in each channel

# Extract node IDs and stream distances
prof_id = list(profiler.data_structure[outlet].values())[0]["ids"] # The IDs of the nodes in the profile. NB only works if main channel only extracted
prof_distances = list(profiler.data_structure[outlet].values())[0]["distances"] # Distances upstream from outlet NB only works if main channel only extracted
prof_xy = np.unravel_index(prof_id, full_shape) # Convert stream nodes into xy coordinates

# Extract downstream predictions along long profile
prof_geochem = sed_norm[prof_id] # Extract geochemistry at each node
# Extract observations along profile
prof_obs = obs_log[np.isin(loc_nodes,prof_id)]
prof_misfit = diff_array[np.isin(loc_nodes,prof_id)]
prof_loc_nodes = loc_nodes[np.isin(loc_nodes,prof_id)] # sample locs on profile
obs_prof_x = np.zeros(prof_obs.size) # Extract stream distance for each locality
for i in np.arange(obs_prof_x.size): # loop sets distance for each locality
    obs_prof_x[i]= prof_distances[prof_id==prof_loc_nodes[i]]
