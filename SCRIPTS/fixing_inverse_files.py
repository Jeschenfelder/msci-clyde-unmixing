import ipywidgets as widgets
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

parser = argparse.ArgumentParser()
parser.add_argument('fieldname', type=str, help='Name of file of faulty inverse')
args = parser.parse_args()
filename = args.fieldname

input_path = '../DATA/INVERSE_RESULTS/' + filename
output_path = '../DATA/INVERSE_RESULTS/fixed_' + filename
field_to_add = np.loadtxt(input_path)

#change nan values to -99
field_to_add = np.nan_to_num(field_to_add,copy=False, nan=-99)

#load in topography:
zr_nc=netCDF4.Dataset('../DATA/Clyde_Topo_100m_working.nc')
zr_ma = zr_nc['z'][:,:]
mg = RasterModelGrid(zr_ma.shape,xy_spacing=(100,100))

zr = mg.add_zeros('node', 'topographic__elevation')
zr += np.reshape(zr_ma.data.astype(float),zr.shape)
dx = 100

# Set up the boundary conditions on the square grid
mg.set_fixed_value_boundaries_at_grid_edges(True,True,True,True)

flat_shape = zr.shape # a tuple to flatten arrays [number of nodes long]
full_shape = mg.shape # the full shape of the grid [rows, columns]

mg.add_field('node', 'inverse_result', field_to_add)
mg.save(output_path, names=['inverse_result'])
