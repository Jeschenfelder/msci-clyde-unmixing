import netCDF4
import sys
import time
import scipy as sp
import landlab
from datetime import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
from landlab import RasterModelGrid
from landlab.components import FlowAccumulator, FastscapeEroder, SinkFillerBarnes, FlowDirectorD8, ChannelProfiler, TrickleDownProfiler
from landlab.components.flow_accum.flow_accum_bw import find_drainage_area_and_discharge
from landlab import imshow_grid
from landlab.utils import get_watershed_mask,get_watershed_outlet,get_watershed_nodes


# Set up grid
print("Setting up model grid...")
zr_nc=netCDF4.Dataset('input_dir/topo_CG.nc')
zr_ma = zr_nc['z'][:,:]
mg = RasterModelGrid(zr_ma.shape,xy_spacing=(200,200))

model_width = mg.shape[1]
model_height = mg.shape[0]

landmask_nc=netCDF4.Dataset('input_dir/landmask.nc')
landmask = landmask_nc['z'][:,:].data.astype(float)

zr = mg.add_zeros('node', 'topographic__elevation')
zr += np.reshape(zr_ma.data.astype(float),zr.shape)
dx = np.loadtxt("input_dir/dx_CG.dat")

# Set up the boundary conditions on the square grid
mg.set_fixed_value_boundaries_at_grid_edges(True,True,True,True)

flat_shape = zr.shape # a tuple to flatten arrays [number of nodes long]
full_shape = mg.shape # the full shape of the grid [rows, columns]

# Sink filling
print("Filling sinks, setting up eroder and flow-routing...")

sfb = SinkFillerBarnes(mg, method='D8',fill_flat=False)
sfb.run_one_step()

# instantiate the flow routing and stream power:

frr = FlowAccumulator(
    mg,
    'topographic__elevation',
    flow_director = 'FlowDirectorD8')
frr.run_one_step()  # flow routing


zrold = zr[mg.nodes] # pre-incision elevations

# instantiate the stream power law
k=3.62
m=0.35
n=1
spr = FastscapeEroder(mg, K_sp=k, m_sp=m, n_sp=n)

dt=0.1# Timestep in Ma
spr.run_one_step(dt)    # erosion: stream power
zrnew = zr[mg.nodes]    # post-incision elevations
incise = zrold - zrnew  # incision per cell
cell_area = mg.dx * mg.dy
qs = incise*cell_area # Volume sediment produced per cell
qsflat = qs.ravel()     # flatten qs for flow routing calculation

a,q = find_drainage_area_and_discharge(mg.at_node['flow__upstream_node_order'], mg.at_node['flow__receiver_node'],runoff = qsflat)
# q is cumulative flux for SPL model
a, q_homo_incis = find_drainage_area_and_discharge(mg.at_node['flow__upstream_node_order'], mg.at_node['flow__receiver_node']) # a is number of nodes
# q is cumulative flux for homogeneous incision

area = mg.at_node['drainage_area']
np.amax(area)
area_threshold = 25 #float(sys.argv[1]) #25 km2
is_drainage = area > (area_threshold*1000000) #km2 to m2
mg.add_field('node','channels',is_drainage,noclobber=False)

# Load in geochemical observations
print("Loading in geochemical observations...")


sample_data = np.loadtxt('input_dir/samples.dat',dtype=str) # [sample #, locality #, x, y]
sample_locs = sample_data[:,2:4].astype(np.float)

channel_xy = np.flip(np.transpose(np.where(is_drainage.reshape(mg.shape))),axis=1)*200 # xy coordinates of channels

nudge = np.zeros(sample_locs.shape) # initiate nudge array

nudge[np.where(sample_data[:,1]=='2'),:] = [-200,-200] # Nudging locality 2 to the SW
nudge[np.where(sample_data[:,1]=='6'),:] = [-600,0] # Nudging locality 6 to the W
nudge[np.where(sample_data[:,1]=='13'),:] = [-200,0] # Nudging locality 13 to the W
nudge[np.where(sample_data[:,1]=='19'),:] = [0,400] # Nudging locality 19 to the N
nudge[np.where(sample_data[:,1]=='20'),:] = [-400,-400] # Nudging locality 20 to the SW
nudge[np.where(sample_data[:,1]=='23'),:] = [-400,0] # Nudging locality 23 to the W
nudge[np.where(sample_data[:,1]=='29'),:] = [-300,0] # Nudging locality 29 to the W
nudge[np.where(sample_data[:,1]=='53'),:] = [-600,-600] # Nudging locality 53 to the SW
nudge[np.where(sample_data[:,1]=='54'),:] = [0,600] # Nudging locality 54 to the N
nudge[np.where(sample_data[:,1]=='56'),:] = [-200,-200] # Nudging locality 56 to the SW
nudge[np.where(sample_data[:,1]=='57'),:] = [-1200,-1200] # Nudging locality 57 to the SW

nudged_locs = sample_locs + nudge # Apply the nudges

# Fit the data to the nearest channel node

fitted_locs = np.zeros(sample_locs.shape) # Initialise snapped locality array
for i in np.arange(nudged_locs.shape[0]):
    sample = nudged_locs[i,:]
    diffs = channel_xy-sample
    distances = np.linalg.norm(diffs,axis=1)
    shortest = np.amin(distances)
    fitted_locs[i,:] = channel_xy[np.where(distances==shortest)]

loc_indxs = np.transpose(np.flip((fitted_locs/200).astype(int),axis=1))
loc_nodes = np.ravel_multi_index(loc_indxs,dims=full_shape)

geochem_raw = np.loadtxt('input_dir/geochem.dat',dtype=str) # Read in data
geochem_raw = np.delete(geochem_raw,[7,53],1) # Delete columns for S and Bi (too many NAs)
elems = geochem_raw[0,1:54] # List of element strings
np.any(geochem_raw=='NA') # Should be False!
obs_data = pd.DataFrame(geochem_raw[1:,],columns=geochem_raw[0,:]) # Cast to DataFrame for quick access
obs_data[elems]=obs_data[elems].astype(float) # Cast numeric data to float

def cost_homo(bdrck_arr,elem):
    """Returns L2norm data misfit for a given bedrock input array (mgkg), calculating predictions using the homogenous incision"""
    a, sed_comp = find_drainage_area_and_discharge(mg.at_node['flow__upstream_node_order'], mg.at_node['flow__receiver_node'],runoff = bdrck_arr)  # composition but homogeneous erosion
    sed_comp_norm = sed_comp/q_homo_incis
    l2norm = np.linalg.norm(np.log10(obs_data[elem]) - np.log10(sed_comp_norm[loc_nodes]))
    return(l2norm)

# Identify active area
print("Identifying coverage area...")

lowest_tay_sample = loc_nodes[sample_data[:,1]=='55'] # Locality 55
tay_catchment = get_watershed_mask(mg,lowest_tay_sample) # extract upstream area of most downstream tay sample

lowest_dee_sample = loc_nodes[sample_data[:,1]=='7'] # Locality 7.
dee_catchment = get_watershed_mask(mg,lowest_dee_sample)

lowest_don_sample = loc_nodes[sample_data[:,1]=='34'] # Locality 34.
don_catchment = get_watershed_mask(mg,lowest_don_sample)

lowest_spey_sample = loc_nodes[sample_data[:,1]=='28'] # Locality 28.
spey_catchment = get_watershed_mask(mg,lowest_spey_sample)

lowest_deveron_sample = loc_nodes[sample_data[:,1]=='40'] # Locality 28.
deveron_catchment = get_watershed_mask(mg,lowest_deveron_sample)

active_area = tay_catchment | dee_catchment | don_catchment | spey_catchment | deveron_catchment

# Defining river outlets to explore geochemical long profiles.

tay_node_eg = mg.grid_coords_to_node_id(100, 400) # a node in tay catchment
tay_outlet = get_watershed_outlet(mg, tay_node_eg) # tay outlet node

dee_node_eg = mg.grid_coords_to_node_id(400, 450) # a node in dee catchment
dee_outlet = get_watershed_outlet(mg, dee_node_eg)

don_node_eg = mg.grid_coords_to_node_id(550, 700) # a node in don catchment
don_outlet = get_watershed_outlet(mg, don_node_eg)

spey_node_eg = mg.grid_coords_to_node_id(500, 450) # a node in spey catchment
spey_outlet = get_watershed_outlet(mg, spey_node_eg)

deveron_node_eg = mg.grid_coords_to_node_id(660, 780) # a node in deveron catchment
deveron_outlet = get_watershed_outlet(mg, deveron_node_eg)

coverage_mask = np.invert(active_area.reshape(full_shape))  # This is a bool array which is used to mask regions with no coverage
coverage_nan_mask = np.zeros(coverage_mask.shape)
coverage_nan_mask[np.where(coverage_mask)] = 'nan'
mg_active_cell_bool = np.where(np.invert(coverage_mask))

#####
# Priors

obs_elems = obs_data[elems]

spey_mth_comp = np.asarray(obs_elems[obs_data['locality'] == '28'])
dee_mth_comp = np.asarray(obs_elems[obs_data['locality'] == '7'])
dev_mth_comp = np.asarray(obs_elems[obs_data['locality'] == '40'])
tay_mth_comp = np.asarray(obs_elems[obs_data['locality'] == '55'])
don_mth_comp = np.asarray(obs_elems[obs_data['locality'] == '34'])

prior_wtd_avg = pd.DataFrame((spey_mth_comp*3012.24 + dee_mth_comp*2009.64 + dev_mth_comp*1407.12 +
           tay_mth_comp*5096.36 + don_mth_comp*1330.4)/(3012.24+2009.64+1407.12+5096.36+1330.4))
prior_wtd_avg.columns = elems
prior_wtd_avg = np.mean(prior_wtd_avg,axis=0)
prior_wtd_avg_log = np.log(prior_wtd_avg)

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

def initiate_blocky_inversion(nx,ny,elem):
    """Initiates an inversion grid for given number of
    cells and element. """
    blox = np.zeros((ny,nx)) + prior_wtd_avg_log[elem] # Set block height to prior in log space
    active_blox = get_active_blocks(nx,ny) # Active cells
    block_x_step = np.ceil(model_width/nx) # Block width
    block_y_step = np.ceil(model_height/ny) # Block height

    return blox,active_blox,block_x_step,block_y_step

# Evaluate results
print("Evaluating results...") 

def abline(slope, intercept):
    # Auxiliary function
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--',c='gray')

def compare_to_average_gbase(input_blocks,input_active_blocks,element,block_x,block_y):
    nx = np.shape(input_blocks)[1]
    ny = np.shape(input_blocks)[0]

    elem_rawdata =  np.loadtxt('DATA/GBASE_MASKED/' + element + '_masked_GBASE.dat',dtype=str)
    good_records = np.invert(elem_rawdata[:,2] == 'NA') # Locations where not NA
    gbase_vals = elem_rawdata[good_records][:,2].astype(float) # log10(elem)
    gbase_locs = elem_rawdata[good_records,:2].astype(float) # x, y
    gbase_locs_inds = np.flip(np.around(gbase_locs/mg.dx).astype(int),axis=1).T # indices. y then x
    gbase_indxs = list(zip(gbase_locs_inds[0,:],gbase_locs_inds[1,:])) # list of tuples of (y,x) indices

    averaged_gbase = np.zeros((ny,nx))
    for i in np.arange(ny):
        for j in np.arange(nx):
            if(active_blocks[i,j]):
                # Boolean array of model cells that correspond to block index (i,j)
                model_index_xrange = (block_width*j,block_width*(j+1))
                model_index_yrange = (block_height*i,block_height*(i+1))
                right_row = np.logical_and(gbase_locs_inds[0,:] >= model_index_yrange[0],gbase_locs_inds[0,:] < model_index_yrange[1])
                right_col = np.logical_and(gbase_locs_inds[1,:] >= model_index_xrange[0],gbase_locs_inds[1,:] < model_index_xrange[1])
                right_row_and_col = np.logical_and(right_row,right_col)
                gbases_in_cell = gbase_locs_inds[:,right_row_and_col]

                if(np.any(gbases_in_cell)):
                    averaged_gbase[i,j] = np.mean(gbase_vals[right_row_and_col])
                else:
                    averaged_gbase[i,j] = 'nan'
            if(not(active_blocks[i,j])):
                averaged_gbase[i,j] = 'nan'

    expanded_pred = expand((np.exp(input_blocks)),block_x,block_y)
    expanded_obs = expand(10**averaged_gbase,block_x,block_y)

    expanded_pred[np.invert(active_area.reshape(full_shape))] = 'nan'
    expanded_obs[np.invert(active_area.reshape(full_shape))] = 'nan'


    plt.subplot(3,2,1)
    plt.imshow(expanded_pred,origin='lower',norm=LogNorm())
    cb = plt.colorbar()
    cb.set_label("Predicted concentration (mg/kg)")
    plt.xlabel("Horizontal Model Cells")
    plt.ylabel("Vertical Model Cells")
    plt.title("Inverse results")

    plt.subplot(3,2,2)
    plt.imshow(expanded_obs,origin='lower',norm=LogNorm())
    cb = plt.colorbar()
    cb.set_label("G-BASE concentration (mg/kg)")
    plt.xlabel("Horizontal Model Cells")
    plt.ylabel("Vertical Model Cells")
    plt.title("G-BASE data")

    pred = np.log10(np.exp(input_blocks[active_blocks]))
    obs = averaged_gbase[active_blocks]
    nan_obs = np.invert(np.isnan(obs)) # Cells with no data
    obs = obs[nan_obs]
    pred = pred[nan_obs]
    misfit = pred-obs
    rms = np.sqrt(np.nanmean(misfit**2))

    plt.subplot(3,2,3)

    plt.plot(np.unique(pred), np.poly1d(np.polyfit(pred, obs, 1))(np.unique(pred)),c='k')
    abline(1,0)
    plt.scatter(x=pred,y=obs,s=4,c=misfit,cmap='seismic',facecolor='grey',norm=TwoSlopeNorm(vcenter=0))
    cb = plt.colorbar()
    cb.set_label("(Underprediction)             Misfit             (Overprediction) ")
    plt.xlabel('Inversion Predictions log10(mg/kg)')
    plt.ylabel(element+'observations log10(mg/kg)')
    ax=plt.gca()
    ax.set_aspect('equal')
    plt.title("Cross-Plot")
    plt.grid(True)

    misfit_spatial = np.log10(np.exp(input_blocks)) - averaged_gbase
    expanded_misfit = expand(misfit_spatial,block_x,block_y)
    expanded_misfit[np.invert(active_area.reshape(full_shape))] = 'nan'
    plt.subplot(3,2,4)
    plt.imshow(expanded_misfit,cmap='seismic',norm=TwoSlopeNorm(vcenter=0),origin='lower')
    cb = plt.colorbar()
    cb.set_label("(Underprediction)             Misfit             (Overprediction) ")
    plt.title("Misfit")
    plt.xlabel("Horizontal Model Cells")
    plt.ylabel("Vertical Model Cells")

    #plt.show()
    plt.tight_layout()

    # Output R2
    correlation_matrix = np.corrcoef(pred, obs)
    correlation_xy = correlation_matrix[0,1]
    r_squared = correlation_xy**2
    print("R-squared value:",r_squared)
    print("RMS misfit =",rms)

    # Writing out files to be read into G-BASE
    expanded_pred = np.log10(expanded_pred)
    expanded_obs = np.log10(expanded_obs)

    expanded_pred[np.invert(active_area.reshape(full_shape))] = -99
    expanded_obs[np.invert(active_area.reshape(full_shape))] = -99
    expanded_misfit[np.invert(active_area.reshape(full_shape))] = -99


    mg.add_field('node','inv_expanded',expanded_pred.reshape(flat_shape),clobber=True) # Expanded model grid
    mg.add_field('node','obs_avg_expanded',expanded_obs.reshape(flat_shape),clobber=True) # Expanded observ grid
    mg.add_field('node','msft_avg_expanded',expanded_misfit.reshape(flat_shape),clobber=True) # Expanded misfit grid
    mg.save('./gmt_tmp/out.asc',names=['inv_expanded','obs_avg_expanded','msft_avg_expanded']) # save as .asc

    np.savetxt("./gmt_tmp/pred_obsavg_misfit.xyz",np.array([pred, obs, misfit]).T) # obs/pred/misfit
    np.savetxt("./gmt_tmp/stats_avg.dat",np.array([rms, r_squared]),fmt='%1.3f') # goodness of fit stats.

def compare_to_gbase(inverse_results_grd,elem):
    """
    Compares the target inverse map to the independent G-BASE dataset.
    Only needs to be provided with an element and a modelled grid and it
    makes the comparison spatially, as a cross plot and shows misfit maps.
    """

    plt.rcParams["figure.figsize"] = (15,15)

    # Set upper and lower bounds
    invmin = np.nanmin(inverse_results_grd+coverage_nan_mask)
    invmax = np.nanmax(inverse_results_grd+coverage_nan_mask)

    # Plot the inverse results
    # plt.subplot(3,2,1)
    # plt.imshow(inverse_results_grd+coverage_nan_mask,origin='lower',norm=LogNorm(vmin=invmin,vmax=invmax))
    # cb = plt.colorbar()
    # cb.set_label("Predicted concentration (mg/kg)")
    # plt.xlabel("Horizontal Model Cells")
    # plt.ylabel("Vertical Model Cells")
    # plt.title("Inverse results")

    # Load in continuous G-BASE map
    tgt_comp_nc=netCDF4.Dataset('input_dir/input_chems/CG_log_'+elem+'.nc')
    tgt_comp_ma = tgt_comp_nc['z'][:,:]
    tgt_comp_bedrock = 10**tgt_comp_ma.data.astype(float).reshape(flat_shape)

    # Plot G-BASE map
    # plt.subplot(3,2,2)
    # plt.imshow(tgt_comp_bedrock.reshape(full_shape)+coverage_nan_mask,origin='lower',norm=LogNorm(vmin=invmin,vmax=invmax))
    # cb = plt.colorbar()
    # cb.set_label("G-BASE concentration (mg/kg)")
    # plt.xlabel("Horizontal Model Cells")
    # plt.ylabel("Vertical Model Cells")
    # plt.title("G-BASE data")

    #### Extract G-BASE point observations and generate cross plot ####
    elem_rawdata =  np.loadtxt('input_dir/masked_'+elem+'.dat',dtype=str)
    good_records = np.invert(elem_rawdata[:,2] == 'NA') # Locations where not NA
    gbase_locs = elem_rawdata[good_records,:2].astype(float) # x, y
    gbase_locs_inds = np.flip(np.around(gbase_locs/mg.dx).astype(int),axis=1).T # indices. y then x

    results_grid = inverse_results_grd.reshape(full_shape).copy()
    predicted = np.log10(results_grid[gbase_locs_inds[0],gbase_locs_inds[1]])
    observed =  elem_rawdata[good_records,2].astype(float)
    misfit = predicted - observed

    # Plot data as a cross-plot
    # plt.subplot(3,2,3)
    # plt.scatter(x=predicted,y=observed,marker='.',s=0.5,color='darkgrey')
    # plt.plot(np.unique(predicted), np.poly1d(np.polyfit(predicted, observed, 1))(np.unique(predicted)),c='k')
    # abline(1,0)
    # plt.xlabel('Inversion Predictions log10(mg/kg)')
    # plt.ylabel(elem+'observations log10(mg/kg)')
    # ax=plt.gca()
    # ax.set_aspect('equal')
    # ax.set_xlim(np.log10([invmin,invmax]))
    # ax.set_ylim(np.log10([invmin,invmax]))
    # plt.grid(True)

    rms_gbase=np.sqrt(np.mean((misfit)**2))
    print("RMS rltv to GBASE =",rms_gbase)
    # Output R2
    correlation_matrix = np.corrcoef(predicted, observed)
    correlation_xy = correlation_matrix[0,1]
    r_squared = correlation_xy**2
    print("R-squared value:", r_squared)


    #### Plotting spatial misfit ####
    # plt.subplot(3,2,4)
    # plt.scatter(x=gbase_locs[:,0],y=gbase_locs[:,1],c=misfit,cmap='seismic',s=4,facecolor='grey',norm=TwoSlopeNorm(vcenter=0))
    # plt.gca().set_aspect('equal')
    # plt.gca().set_facecolor('grey')
    # cb = plt.colorbar()
    # cb.set_label("(Underprediction)                  Misfit                  (Overprediction) ")

    # ### Plotting distribution of misfits ####
    # plt.subplot(3,2,5)
    # plt.hist(x=misfit,bins=100)
    # plt.xlabel('(Underprediction)                  Misfit                  (Overprediction) ')
    # plt.tight_layout()
    #plt.show()

    # Writing out files to be read into G-BASE
    results_grid = np.log10(results_grid).reshape(full_shape)
    tgt_comp_bedrock = np.log10(tgt_comp_bedrock).reshape(full_shape)

    results_grid[np.invert(active_area.reshape(full_shape))] = -99
    tgt_comp_bedrock[np.invert(active_area.reshape(full_shape))] = -99

    mg.add_field('node','inv_expanded',results_grid.reshape(flat_shape),clobber=True) # Expanded model grid
    mg.save('./gmt_tmp/out_fullres_inversion.asc',names=['inv_expanded']) # save as .asc


    np.savetxt("./gmt_tmp/full_res_spatialmisfit.xyz",np.array([gbase_locs[:,0],gbase_locs[:,1],misfit]).T) # x/y/misfit
    np.savetxt("./gmt_tmp/pred_obsfullres_misfit.xyz",np.array([predicted, observed, misfit]).T) # obs/pred/misfit
    np.savetxt("./gmt_tmp/stats_fullres.dat",np.array([rms_gbase, r_squared]),fmt='%1.3f') # goodness of fit stats.

def evaluate_downstream_samps(results_grid,elem):
    """ Produce a cross plot for the predicted and observed
    downstream geochemical samples"""

    plt.rcParams["figure.figsize"] = (5,5)
    comp_bedrock = results_grid.reshape(flat_shape)

    # For homogeneous incision
    a, sed_comp = find_drainage_area_and_discharge(mg.at_node['flow__upstream_node_order'], mg.at_node['flow__receiver_node'],runoff = comp_bedrock)  # composition but homogeneous erosion
    with np.errstate(divide='ignore'):
        with np.errstate(invalid='ignore'):
            sed_comp_norm = sed_comp /q_homo_incis
    sed_comp_norm[q_homo_incis==0] = comp_bedrock[q_homo_incis==0]
    preds = sed_comp_norm[loc_nodes]
    obs = obs_data[elem]

    misfit = np.log10(preds)-np.log10(obs)

    rms=np.sqrt(np.mean((misfit)**2))
    print("RMS downstream = ",rms)
    correlation_matrix = np.corrcoef(preds, obs)
    correlation_xy = correlation_matrix[0,1]
    r_squared = correlation_xy**2
    print("R-squared value downstream:", r_squared)

    plt.scatter(x=preds,y= obs, c=misfit, s=40,cmap='seismic',norm=TwoSlopeNorm(vcenter=0))
    cb = plt.colorbar()
    cb.set_label("(Underprediction)                  Misfit                  (Overprediction) ")
    # Add regression line
    plt.plot(np.unique(preds), 10**np.poly1d(np.polyfit(np.log10(preds), np.log10(obs), 1))(np.unique(np.log10(preds))),c='k')
    plt.xlabel('Predicted downstream '+elem+', mg/kg')
    plt.ylabel('Observed downstream'+elem+', mg/kg')
    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.axis('equal')
    abline(1,0)
    #plt.show()

    np.savetxt("./gmt_tmp/downstream_pred_obs_misfit.xyz",np.array([np.log10(preds), np.log10(obs), misfit]).T) # obs/pred/misfit
    np.savetxt("./gmt_tmp/stats_downstream.dat",np.array([rms, r_squared]),fmt='%1.3f') # goodness of fit stats.

def misfit_map_downstream(results_grid,elem):
    """ Produce a spatial map of the misfit
    for the downstream samples"""

    plt.rcParams["figure.figsize"] = (7,5)
    comp_bedrock = results_grid.reshape(flat_shape)

    # For homogeneous incision
    a, sed_comp = find_drainage_area_and_discharge(mg.at_node['flow__upstream_node_order'], mg.at_node['flow__receiver_node'],runoff = comp_bedrock)  # composition but homogeneous erosion
    with np.errstate(divide='ignore'):
        with np.errstate(invalid='ignore'):
            sed_comp_norm = sed_comp /q_homo_incis
    sed_comp_norm[q_homo_incis==0] = comp_bedrock[q_homo_incis==0]
    preds = sed_comp_norm[loc_nodes]
    obs = obs_data[elem]

    misfit = np.log10(preds)-np.log10(obs)

    # greys_cmap = plt.cm.Greys
    # greys_cmap.set_bad(color='darkgrey')
    # plt.imshow(is_drainage.reshape(full_shape)+landmask,origin='lower',cmap='Greys')
    # plt.scatter(x=fitted_locs[:,0]/mg.dx, y=fitted_locs[:,1]/mg.dx, c=misfit, s=60,cmap='seismic',norm=TwoSlopeNorm(vcenter=0))
    # cb = plt.colorbar()
    # cb.set_label("(Underprediction)                  Misfit                  (Overprediction) ")
    # plt.tight_layout()
    # #plt.show()

    nanned_channels = is_drainage.reshape(full_shape)
    nanned_channels[np.invert(active_area.reshape(full_shape))] = -99
    mg.add_field('node','nanned_channels',nanned_channels,noclobber=False)
    mg.add_field('node','outsedcomp',np.log10(sed_comp_norm),noclobber=False)

    mg.save('./gmt_tmp/out_nanned_channels.asc',names=['nanned_channels']) # save channels as .asc
    mg.save('./gmt_tmp/out_sedcomp.asc',names=['outsedcomp']) # save channels as .asc

    np.savetxt("./gmt_tmp/downstream_spatialmisfit.xyz",np.array([fitted_locs[:,0], fitted_locs[:,1], misfit]).T) # obs/pred/misfit

def compare_long_profiles(results_grid,elem):
    """ Produce long profiles for each major catchment's geochemistry
     overlain with the observations"""

    plt.rcParams["figure.figsize"] = (15,5)

    comp_bedrock = results_grid.reshape(flat_shape)
    # Calculate sediment composition
    a, sed_comp = find_drainage_area_and_discharge(mg.at_node['flow__upstream_node_order'], mg.at_node['flow__receiver_node'],runoff = comp_bedrock)  # composition but homogeneous erosion
    with np.errstate(divide='ignore'):
        with np.errstate(invalid='ignore'):
            sed_comp_norm = sed_comp /q_homo_incis
    sed_comp_norm[q_homo_incis==0] = comp_bedrock[q_homo_incis==0]
    preds = sed_comp_norm[loc_nodes] # Extract predictions
    obs = obs_data[elem] # Extract observations
    misfit = np.log10(preds)-np.log10(obs)
    misfitmin,misfitmax = np.amin(misfit),np.amax(misfit)

    plt.figure()
    ax1 = plt.subplot(1,2,1) # Long profile figure
    ax2 = plt.subplot(1,2,2) # Map

    riverlist = ['Spey','Tay','Dee','Don','Deveron']
    for river in riverlist:
        if(river=='Spey'):
            outlet=spey_outlet
        if(river=='Tay'):
            outlet=tay_outlet
        if(river=='Dee'):
            outlet=dee_outlet
        if(river=='Don'):
            outlet=don_outlet
        if(river=='Deveron'):
            outlet=deveron_outlet

        profiler = ChannelProfiler(mg,main_channel_only=True,outlet_nodes=[outlet])
        profiler.run_one_step() # Extract longest profile in each channel

        # Extract node IDs and stream distances
        prof_id = list(profiler.data_structure[outlet].values())[0]["ids"] # NB only works if main channel only extracted
        prof_distances = list(profiler.data_structure[outlet].values())[0]["distances"] # NB only works if main channel only extracted
        prof_xy = np.unravel_index(prof_id, full_shape) # Convert stream nodes into xy coordinates

        # Extract downstream predictions along long profile
        prof_geochem = sed_comp_norm[prof_id]
        # Extract observations along profile
        obs_prof_val = obs[np.isin(loc_nodes,prof_id)]
        obs_prof_misfit = misfit[np.isin(loc_nodes,prof_id)] # Calculate misfit
        prof_loc_nodes = loc_nodes[np.isin(loc_nodes,prof_id)] # sample locs on profile
        obs_prof_x = np.zeros(obs_prof_val.size) # Extract stream distance for each locality
        for i in np.arange(obs_prof_x.size): # loop is used to deal with duplicates
            obs_prof_x[i]= prof_distances[prof_id==prof_loc_nodes[i]]

        # Plot geochemical profile
        ax1.plot(prof_distances/1e3,prof_geochem)
        ax1.scatter(obs_prof_x/1e3,obs_prof_val)#,c=obs_prof_misfit,s=60,cmap='seismic',vmin = misfitmin,vmax=misfitmax, norm=TwoSlopeNorm(vcenter=0))

        # Plot map of long profiles
        ax2.imshow(zr.reshape(full_shape)+landmask,cmap='terrain',origin='lower',norm=TwoSlopeNorm(vcenter=0))
        ax2.scatter(prof_xy[1],prof_xy[0],s=2)
        ax2.scatter(x=fitted_locs[np.isin(loc_nodes,prof_id),0]/mg.dx,
                y=fitted_locs[np.isin(loc_nodes,prof_id),1]/mg.dx,s=30,marker='x',c='black')

        np.savetxt("./gmt_tmp/"+river+".xy",np.array([prof_xy[1]*mg.dx,prof_xy[0]*mg.dx]).T) #
        np.savetxt("./gmt_tmp/"+river+"_distances_predconc.xy",np.array([prof_distances/1e3,np.log10(prof_geochem)]).T) #
        np.savetxt("./gmt_tmp/"+river+"_distances_obsconc.xy",np.array([obs_prof_x/1e3,np.log10(obs_prof_val)]).T) #


    # Tidy up figures
    ax1.set_xlabel('Distance to mouth / km')
    ax1.set_yscale('log')
    ax1.set_ylabel('Predicted downstream'+elem+', mg/kg')
    ax1.legend(labels=riverlist)
    ax2.set_xlabel('Horizontal Model Cells')
    ax2.set_ylabel('Vertical Model Cells')

    plt.tight_layout()
    #plt.show()

# Perform the comparison
element  = sys.argv[3] #'U'
nx = int(sys.argv[4]) #40
ny = int(sys.argv[5]) #34
infile_interp = sys.argv[1] # 'HPC_inverse_results/20x20_U_nonlin_-1.0_interp.npy'
infile_x = sys.argv[2] # 'HPC_inverse_results/20x20_U_nonlin_-1.0_interp.npy'

input_grd = np.load(infile_interp)
input_x = np.load(infile_x)

blocks,active_blocks,block_width,block_height = initiate_blocky_inversion(nx,ny,element)
blocks[active_blocks] = input_x

compare_to_average_gbase(blocks,active_blocks,element,block_width,block_height)
compare_to_gbase(input_grd,element)
evaluate_downstream_samps(input_grd,element)
misfit_map_downstream(input_grd,element)
compare_long_profiles(input_grd,element)
