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

########################################### setting up the topographic data ######################################
zr_nc=netCDF4.Dataset('../DATA/Clyde_Topo_100m_working.nc')
zr_ma = zr_nc['z'][:,:]
#save array as filled_topo.npy
topography = np.load('filled_topography.npy')
#add array as topography field
mg = RasterModelGrid(zr_ma.shape,xy_spacing=(100,100))
zr = mg.add_field('node', 'topographic__elevation', topography)
dx = 100

flat_shape = zr.shape # a tuple to flatten arrays [number of nodes long]
full_shape = mg.shape # the full shape of the grid [rows, columns]

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
sample_data = np.loadtxt('../DATA/filtered_sample_loc.dat',dtype=str) # [x, y, sample #]
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

obs_data = pd.read_csv('../DATA/converted_chem_data.csv') #read in geochem data
elems =  obs_data.columns[1:].tolist() # List of element strings
obs_data[elems]=obs_data[elems].astype(float) # Cast numeric data to float

#loading in compositional data:
obs_elems = obs_data[elems]
clyde_mth_comp = np.asarray(obs_elems[obs_data['SAMPLE_No'] == 700012])

#average 
prior_wtd_avg = pd.DataFrame(clyde_mth_comp)
prior_wtd_avg.columns = elems

prior_wtd_avg = np.mean(prior_wtd_avg,axis=0)
prior_wtd_avg_log = np.log(prior_wtd_avg)

#load in upstream area for each locality:
unique_locs = np.unique(sample_data[:,2])
loc_areas = np.load('loc_areas.npy')

#set dimensions and lowest sample
model_width = mg.shape[1] # number of x cells in topographic model grid
model_height = mg.shape[0] # number of y cells in topographic model grid
lowest_sample = loc_nodes[sample_data[:,2]== '700012'] # Locality 700012

#load in active blocks:
active_blocks = np.load('active_blocks_84x74.npy')
active_area = get_watershed_mask(mg,lowest_sample) # extract upstream area of most downstream Clyde sample 

#set up block dimensions:
nx=84 #  # <<<<<<<<<<<<<<<<   Change number of x blocks in inversion grid
ny=74 #  # <<<<<<<<<<<<<<<<   Change number of y blocks in inversion grid
block_width = np.ceil(model_width/nx) # set block width
block_height = np.ceil(model_height/ny) # set block height

############################# Load in masked G-BASE and average out by block ######################

#define expand function:
def expand(block_grid,block_x,block_y):
    """Expands low res array of block heights into 
    model grid array that can be fed into topographic
    model. Note that blocks at the upper and eastern 
    perimeter are clipped if number of blocks doesn't 
    divide number of model cells. 
    
    block_x and block_y are the number of model cells 
    in each block in x and y dir respectively"""
    return(block_grid.repeat(block_y, axis=0).repeat(block_x, axis=1)[:model_height,:model_width])

def compare_to_average_gbase(input_blocks,input_active_blocks,element,block_x = block_width,block_y = block_height, num_x = nx,num_y = ny):

    elem_rawdata =  np.loadtxt('DATA/GBASE_MASKED/' + element + '_masked_GBASE.dat',dtype=str) #loading in G-BASE points inside active area
    #good_records = np.invert(elem_rawdata[:,2] == 'NA') # Locations where not NA
    gbase_vals = np.log10(elem_rawdata[:,2].astype(float)) # log10(elem)
    gbase_locs = elem_rawdata[:,:2].astype(float) # x, y
    gbase_locs_inds = np.flip(np.around(gbase_locs/mg.dx).astype(int),axis=1).T # indices. y then x
    gbase_indxs = list(zip(gbase_locs_inds[0,:],gbase_locs_inds[1,:])) # list of tuples of (y,x) indices

    averaged_gbase = np.zeros((num_y,num_x))
    for i in np.arange(num_y):
        for j in np.arange(num_x):
            if(active_blocks[i,j]):
                # Boolean array of model cells that correspond to block index (i,j)
                model_index_xrange = (block_width*j,block_width*(j+1))
                model_index_yrange = (block_height*i,block_height*(i+1))
                right_row = np.logical_and(gbase_locs_inds[0,:] >= model_index_yrange[0],gbase_locs_inds[0,:] < model_index_yrange[1])
                right_col = np.logical_and(gbase_locs_inds[1,:] >= model_index_xrange[0],gbase_locs_inds[1,:] < model_index_xrange[1])
                right_row_and_col = np.logical_and(right_row,right_col)
                gbases_in_cell = gbase_locs_inds[:,right_row_and_col]

                if(np.any(gbases_in_cell)):
                    averaged_gbase[i,j] = np.nanmean(gbase_vals[right_row_and_col])
                else:
                    averaged_gbase[i,j] = 'nan'
            if(not(active_blocks[i,j])):
                averaged_gbase[i,j] = 'nan'

    #expand averaged G-BASE to full grid and add NAN points to outside of active area
    expanded_obs = expand(10**averaged_gbase,block_x,block_y)
    expanded_obs[np.invert(active_area.reshape(full_shape))] = -99

    #save expanded average as .asc file:
    output_path = 
    field_to_add = mg.add_field('nodes', 'average_gbase', expanded_obs)
    