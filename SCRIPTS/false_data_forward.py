import ipywidgets as widgets
import netCDF4
import sys
import time
from numpy.core.fromnumeric import amin
from numpy.core.numeric import full
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
tension=0.25 #<<<<<<<<<<<<<<<<<<<<<<<< change to correct tension factor
interpolate_input = 'DATA/INTERPOLATED_GBASE/gbase_log_' + element + '_T' + str(tension) + '.nc' #path to interpolated G-BASE data
cairngorms_interpolate_input = 'DATA/INTERPOLATED_GBASE/cairngorms_gbase_log_' + element + '_T' + str(tension) + '.nc'
rand_result_output_path = 'DATA/FORWARDMODEL_RESULTS/' + element + '_random_gbase_log_sed.asc' #path to full saved output
cairngorms_result_output_path = 'DATA/FORWARDMODEL_RESULTS/' + element + '_cairngorms_gbase_log_sed.asc' #path to full saved output
flipped_results_output_path = 'DATA/FORWARDMODEL_RESULTS/' + element + '_flipped_gbase_log_sed.asc' #path to full saved output
rand_misfit_output_path = 'DATA/FORWARDMODEL_RESULTS/' + element + '_obs_v_pred_random.txt' #path to output at observed localities
cairngorms_misfit_output_path = 'DATA/FORWARDMODEL_RESULTS/' + element + '_obs_v_pred_cairngorms.txt' #path to output at observed localities
flipped_misfit_output_path = 'DATA/FORWARDMODEL_RESULTS/' + element + '_obs_v_pred_flipped.txt' #path to output at observed localities
path_obs_profile = 'DATA/FORWARDMODEL_RESULTS/' + element +'_obs_profile.txt'
rand_path_pred_profile = 'DATA/FORWARDMODEL_RESULTS/' + element +'_pred_profile_random.txt'
flipped_path_pred_profile = 'DATA/FORWARDMODEL_RESULTS/' + element +'_pred_profile_flipped.txt'
cairngorms_path_pred_profile = 'DATA/FORWARDMODEL_RESULTS/' + element +'_pred_profile_cairngorms.txt'
rand_path_gbase = 'DATA/INTERPOLATED_GBASE/random_gbase_log_' + element + '_T0.25.asc' #path to save the randomisesd G-BASE data
flipped_path_gbase = 'DATA/INTERPOLATED_GBASE/flipped_gbase_log_' + element + '_T0.25.asc' #path to save the flipped G-BASE data
cairngorms_path_gbase = 'DATA/INTERPOLATED_GBASE/cairngorms_gbase_log_' + element + '_T0.25.asc' #path to save the randomisesd G-BASE data

#loading in filled topographic data:
zr_nc=netCDF4.Dataset('DATA/Clyde_Topo_100m_working.nc')
zr_ma = zr_nc['z'][:,:]
topography = np.load('SCRIPTS/filled_topography.npy')

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

#Find total sediment flux first, assuming homogenous incision:
a, q = find_drainage_area_and_discharge(mg.at_node['flow__upstream_node_order'], mg.at_node['flow__receiver_node']) # a is number of nodes

#load in and create flipped G-BASE and run model:
comp_nc=netCDF4.Dataset(interpolate_input)
comp_ma = comp_nc['z'][:,:]
comp_ma = comp_ma.reshape(flat_shape) #change to 1D array prior to shuffling
comp_flipped = np.flip(comp_ma)
comp_flipped = comp_flipped.reshape(full_shape) #reform 2D array
comp_flip = mg.add_zeros('node','flip_bdrck')
comp_flip += np.reshape(((10**comp_flipped).data.astype(float)),comp_flip.shape) # Convert to raw values from log10
comp_flip_log = mg.add_zeros('node','log_flip_bdrck')
comp_flip_log += np.reshape((comp_flipped.data.astype(float)),comp_flip_log.shape)
#mg.save(flipped_path_gbase, names=['log_flip_bdrck']) #save the randomised G-BASE data

#randomise G-BASE data:
np.random.shuffle(comp_ma) #shuffle array to randomise
comp_ma = comp_ma.reshape(full_shape) #reform 2D array
comp_rand_log = mg.add_zeros('node','log_rand_bdrck')
comp_rand_log += np.reshape((comp_ma.data.astype(float)),comp_rand_log.shape)
comp_rand = mg.add_zeros('node','rand_bdrck')
comp_rand += np.reshape(((10**comp_ma).data.astype(float)),comp_rand.shape) # Convert to raw values from log10
#mg.save(rand_path_gbase, names=['log_rand_bdrck']) #save the randomised G-BASE data

#load in cairngorms G-BASE data:
comp_nc=netCDF4.Dataset(cairngorms_interpolate_input)
comp_ma = comp_nc['z'][:,:]
comp_ma = comp_ma[:full_shape[0], :full_shape[1]]
comp_cairngorms_log = mg.add_zeros('node','log_cairngorms_bdrck')
comp_cairngorms_log += np.reshape((comp_ma.data.astype(float)),comp_cairngorms_log.shape)
comp_cairngorms = mg.add_zeros('node','cairngorms_bdrck')
comp_cairngorms += np.reshape(((10**comp_ma).data.astype(float)),comp_rand.shape) # Convert to raw values from log10
mg.save(cairngorms_path_gbase, names=['log_cairngorms_bdrck'])


#Run forward model using composition and homogenous erosion for random data:
a, sed_comp = find_drainage_area_and_discharge(mg.at_node['flow__upstream_node_order'], mg.at_node['flow__receiver_node'],runoff = comp_rand)
sed_norm_rand = sed_comp/q #normalise composition by total sediment flux
sed_norm_rand[q==0] = comp_rand[q==0] #setting composition to bedrock composition where sed flux is 0
#visualise by turning back to log10 and running through channel system:
sed_comp_norm_rand_channel = np.log10(sed_norm_rand) * is_drainage

#Add model results to raster model grid:
mg.add_field('node','homo_incis_sed_rand',sed_norm_rand,noclobber=False)
mg.add_field('node','homo_incis_log_sed_rand',np.log10(sed_norm_rand),noclobber=False)
mg.add_field('node','homo_incis_log_sed_channel_rand',sed_comp_norm_rand_channel,noclobber=False)

#saving the result:
#mg.save(rand_result_output_path, names=['homo_incis_log_sed_channel_rand'])

#Run forward model using composition and homogenous erosion for cairngorms data:
a, sed_comp = find_drainage_area_and_discharge(mg.at_node['flow__upstream_node_order'], mg.at_node['flow__receiver_node'],runoff = comp_cairngorms)
sed_norm_cairngorms = sed_comp/q #normalise composition by total sediment flux
sed_norm_cairngorms[q==0] = comp_cairngorms[q==0] #setting composition to bedrock composition where sed flux is 0
#visualise by turning back to log10 and running through channel system:
sed_comp_norm_cairngorms_channel = np.log10(sed_norm_cairngorms) * is_drainage

#Add model results to raster model grid:
mg.add_field('node','homo_incis_sed_cairngorms',sed_norm_cairngorms,noclobber=False)
mg.add_field('node','homo_incis_log_sed_cairngorms',np.log10(sed_norm_cairngorms),noclobber=False)
mg.add_field('node','homo_incis_log_sed_channel_cairngorms',sed_comp_norm_cairngorms_channel,noclobber=False)

#saving the result:
#mg.save(cairngorms_result_output_path, names=['homo_incis_log_sed_channel_cairngorms'])

#Run forward model using composition and homogenous erosion for cairngorms data:
a, sed_comp = find_drainage_area_and_discharge(mg.at_node['flow__upstream_node_order'], mg.at_node['flow__receiver_node'],runoff = comp_flip)
sed_norm_flip = sed_comp/q #normalise composition by total sediment flux
sed_norm_flip[q==0] = comp_flip[q==0] #setting composition to bedrock composition where sed flux is 0
#visualise by turning back to log10 and running through channel system:
sed_comp_norm_flip_channel = np.log10(sed_norm_flip) * is_drainage

#Add model results to raster model grid:
mg.add_field('node','homo_incis_sed_flip',sed_norm_flip,noclobber=False)
mg.add_field('node','homo_incis_log_sed_flip',np.log10(sed_norm_flip),noclobber=False)
mg.add_field('node','homo_incis_log_sed_channel_flip',sed_comp_norm_flip_channel,noclobber=False)

#saving the result:
#mg.save(flipped_results_output_path, names=['homo_incis_log_sed_channel_flip'])

####################################calculating data misfit between observations and predictions:######################
sample_data = np.loadtxt('DATA/filtered_sample_loc_bad_points.dat',dtype=str) # [x, y, sample #]
sample_locs = sample_data[:,0:2].astype(float)
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

obs_data = pd.read_csv('DATA/converted_chem_data_bad_points.csv') #read in geochem data
elems =  obs_data.columns[1:].tolist() # List of element strings
obs_data[elems]=obs_data[elems].astype(float) # Cast numeric data to float

#creating arrays of observations and predictions:
obs_log = np.log10(obs_data[element])
rand_pred_log = sed_comp_norm_rand_channel[loc_nodes] #get random point observations
cairngorms_pred_log = sed_comp_norm_cairngorms_channel[loc_nodes]
flip_pred_log = sed_comp_norm_flip_channel[loc_nodes]
rand_diff_array = obs_log-rand_pred_log
cairngorms_diff_array = obs_log-cairngorms_pred_log
flip_diff_array = obs_log -flip_pred_log
rand_out_array = np.array([sample_data[:,2].astype(int), sample_data[:,0].astype(float), sample_data[:,1].astype(float), rand_pred_log, obs_log, rand_diff_array]).T #create output array that has the form sample#, x, y, obs, pred, obs-pred
cairngorms_out_array = np.array([sample_data[:,2].astype(int), sample_data[:,0].astype(float), sample_data[:,1].astype(float), cairngorms_pred_log, obs_log, cairngorms_diff_array]).T #create output array that has the form sample#, x, y, obs, pred, obs-pred
flip_out_array = np.array([sample_data[:,2].astype(int), sample_data[:,0].astype(float), sample_data[:,1].astype(float), flip_pred_log, obs_log, flip_diff_array]).T #create output array that has the form sample#, x, y, obs, pred, obs-pred
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
rand_prof_geochem = sed_norm_rand[prof_id] # Extract geochemistry at each node for random data
cairngorms_prof_geochem = sed_norm_cairngorms[prof_id]
flip_prof_geochem = sed_norm_flip[prof_id]
# Extract observations along profile
prof_obs = obs_log[np.isin(loc_nodes,prof_id)]
rand_prof_misfit = rand_diff_array[np.isin(loc_nodes,prof_id)] #get misfit for random points
cairngorms_prof_misfit = cairngorms_diff_array[np.isin(loc_nodes,prof_id)] #get misfit for cairngorms points
flip_prof_misfit = flip_diff_array[np.isin(loc_nodes,prof_id)] #get misfit for cairngorms points
prof_loc_nodes = loc_nodes[np.isin(loc_nodes,prof_id)] # sample locs on profile
obs_prof_x = np.zeros(prof_obs.size) # Extract stream distance for each locality
for i in np.arange(obs_prof_x.size): # loop sets distance for each locality
    obs_prof_x[i]= prof_distances[prof_id==prof_loc_nodes[i]]

#saving all profile distances and chemistry
obs_profile_output = np.array([obs_prof_x,prof_obs]).T
rand_pred_profile_output = np.array([prof_distances,np.log10(rand_prof_geochem)]).T #converting random geochem into log10
cairngorms_pred_profile_output = np.array([prof_distances,np.log10(cairngorms_prof_geochem)]).T #converting cairngorms geochem into log10
flip_pred_profile_output = np.array([prof_distances,np.log10(flip_prof_geochem)]).T #converting cairngorms geochem into log10
np.savetxt(path_obs_profile, obs_profile_output) #save observations
#save random data:
np.savetxt(rand_path_pred_profile, rand_pred_profile_output) #
np.savetxt(rand_misfit_output_path, rand_out_array, fmt = ['%d', '%.18f', '%.18f','%.5f', '%.5f', '%.5f'], header='SAMPLE_No X_COORD Y_COORD CONC_PREDICTION CONC_OBSERVATION MISFIT')
#save cairngorms data:
np.savetxt(cairngorms_path_pred_profile, cairngorms_pred_profile_output)
np.savetxt(cairngorms_misfit_output_path, cairngorms_out_array, fmt = ['%d', '%.18f', '%.18f','%.5f', '%.5f', '%.5f'], header='SAMPLE_No X_COORD Y_COORD CONC_PREDICTION CONC_OBSERVATION MISFIT')
#save flipped data:
np.savetxt(flipped_path_pred_profile, flip_pred_profile_output)
np.savetxt(flipped_misfit_output_path, flip_out_array, fmt = ['%d', '%.18f', '%.18f','%.5f', '%.5f', '%.5f'], header='SAMPLE_No X_COORD Y_COORD CONC_PREDICTION CONC_OBSERVATION MISFIT')


'''
# Plot geochemical profile
fig, (ax1, ax2) = plt.subplots(2,1)
ax1.plot(prof_distances/1e3,np.log10(prof_geochem))
ax1.scatter(obs_prof_x/1e3,prof_obs)#,c=obs_prof_misfit,s=60,cmap='seismic',vmin = misfitmin,vmax=misfitmax, norm=TwoSlopeNorm(vcenter=0))

# Plot map of long profiles
ax2.imshow(zr.reshape(full_shape),cmap='terrain',origin='lower')#,norm=TwoSlopeNorm(vcenter=0))
ax2.scatter(prof_xy[1],prof_xy[0],s=2)
ax2.scatter(x=fitted_locs[np.isin(loc_nodes,prof_id),0]/mg.dx,
        y=fitted_locs[np.isin(loc_nodes,prof_id),1]/mg.dx,s=30,marker='x',c='black')
plt.show()
'''