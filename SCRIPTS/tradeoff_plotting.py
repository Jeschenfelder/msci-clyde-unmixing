import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import pandas as pd

elem = 'Mg' #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< change to correct element!
infile = 'DATA/INVERSE_RESULTS/' + elem + '_results/' + elem + '_all_roughness_misfit.txt'

df = pd.read_csv(infile,delimiter=' ', header=None, dtype=float) # read in data
lambdas = df.iloc[0].to_list() # store all lambda values
misfit = df.iloc[3].to_list() # store all misfits
roughness = (df.iloc[1] + df.iloc[2]).to_list() # calculating total roughness

#create lambda colourmap:


#plot trade off curve:
plt.scatter(roughness, misfit, c = )
plt.show()
