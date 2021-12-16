import ipywidgets as widgets
import netCDF4
import sys
import time
from numpy.core.numeric import isclose
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
parser.add_argument('element', type=str, help='Characteristic of synthetic')
parser.add_argument('fieldname', type=str, help='Name of inverse .npy file')
args = parser.parse_args()
elem  = args.element
filename = args.fieldname

# define input and output paths:
input_path = '../../DATA/SYNTH_RESULTS/' + elem +'_results/' + filename
output_path = '../../DATA/SYNTH_RESULTS/' + elem +'_results/' + elem + '_inverse_output.asc'

#loading in inverse results, active area and setting inactives to NaN:
field_to_add = np.load(input_path).astype(float)
print(np.shape(field_to_add))


#load in topography:
zr_nc=netCDF4.Dataset('../../DATA/Clyde_Topo_100m_working.nc')
zr_ma = zr_nc['z'][:,:]
mg = RasterModelGrid(zr_ma.shape,xy_spacing=(100,100))

zr = mg.add_zeros('node', 'topographic__elevation')
zr += np.reshape(zr_ma.data.astype(float),zr.shape)
dx = 100

# Set up the boundary conditions on the square grid
mg.set_fixed_value_boundaries_at_grid_edges(True,True,True,True)

#Set dimensions for output:
flat_shape = zr.shape # a tuple to flatten arrays [number of nodes long]
full_shape = mg.shape # the full shape of the grid [rows, columns]
model_width = mg.shape[1] # number of x cells in topographic model grid
model_height = mg.shape[0] # number of y cells in topographic model grid

block_width = np.ceil(model_width/np.shape(field_to_add)[1]) # calculate final block width from model width and number of columns in the inverse model
block_height = np.ceil(model_height/np.shape(field_to_add)[0]) # calculate final block height from model height and number of rows in the inverse model


field_to_add = np.nan_to_num(field_to_add,copy=False, nan=-99)

mg.add_field('node', 'inverse_result', field_to_add)
mg.save(output_path, names=['inverse_result'])
