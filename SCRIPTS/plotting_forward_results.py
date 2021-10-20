import numpy as np
from typing import no_type_check
from landlab.components import ChannelProfiler
from landlab import RasterModelGrid
from landlab.components import FlowAccumulator, FastscapeEroder, SinkFillerBarnes, FlowDirectorD8, ChannelProfiler, TrickleDownProfiler
from landlab.components.flow_accum.flow_accum_bw import find_drainage_area_and_discharge
from landlab import imshow_grid
from landlab.utils import get_watershed_mask,get_watershed_outlet,get_watershed_nodes
import netCDF4

#inputs:
element='Mg' #<<<<<<<<<<<<<<<<<<<<<<<< change to correct element
results_path = 'msci-clyde-unmixing/DATA/FORWARDMODEL_RESULTS/' + element + '_obs_v_pred.txt'

# load in topography:
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

#load in results from model
result_array = np.loadtxt(results_path)
