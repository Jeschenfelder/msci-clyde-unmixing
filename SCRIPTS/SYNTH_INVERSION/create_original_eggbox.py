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

eggbox_size_m = int(sys.argv[1]) #enter size of eggbox
output_path = '../../DATA/SYNTH_RESULTS/' + str(int(eggbox_size_m/1000)) + 'km_results/original_eggbox.asc' #specifying output path
#### Set up the model grid with input topography ####
zr_nc=netCDF4.Dataset('../../DATA/Clyde_Topo_100m_working.nc')
zr_ma = zr_nc['z'][:,:]
mg = RasterModelGrid(zr_ma.shape,xy_spacing=(100,100))

zr = mg.add_zeros('node', 'topographic__elevation')
zr += np.reshape(zr_ma.data.astype(float),zr.shape)
dx = 100
eggbox_size=float(eggbox_size_m)/dx

# Set up the boundary conditions on the square grid
mg.set_fixed_value_boundaries_at_grid_edges(True,True,True,True)

flat_shape = zr.shape # a tuple to flatten arrays [number of nodes long]
full_shape = mg.shape # the full shape of the grid [rows, columns]
model_width = mg.shape[1] # number of x cells in topographic model grid
model_height = mg.shape[0] # number of y cells in topographic model grid

def eggbox(w,max=1e4,min=1e1):
    """Returns a sin wave chequerboard"""
    out = np.zeros(full_shape)
    for i in np.arange(out.shape[0]):
        for j in np.arange(out.shape[1]):
            out[i,j] += ((((np.sin(i*np.pi/w) * np.sin(j*np.pi/w))+1)/2)*(max-min) +min)
    return(out)

# Generate downstream data and priors
source_region = eggbox(eggbox_size)

plt.imshow(source_region, origin='lower')
plt.show()

mg.add_field('node', 'original_eggbox', source_region)
mg.save(output_path, names=['original_eggbox'])
