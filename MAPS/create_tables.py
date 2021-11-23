import pandas as pd
from pyproj import Proj, transform

#prepping tranformation:
proj_in = Proj('epsg:4326')
proj_out = Proj('epsg:27700')

df = pd.read_csv('DATA/filtered_sample_loc.dat', delimiter='\t', header=None) #load in filtered locs as x, y, Sample#
df_upper = pd.read_csv('DATA/ClydeHighOrderDrainage.dat', delimiter= '\t') #load in upper clyde data
df_estuary = pd.read_csv('DATA/Compiled_Clyde_Surface.dat') #load in estuarine data

#creating output dataframe from upper:
df_final = df_upper[['SAMPLE_No', 'Easting', 'Northing']].copy() #all upper samples used
df_final = df_final.drop(0)
print(df_final)

#filter estuary data:
needed_samples = df[2].to_list() #get all used sample#
df_filtered = df_estuary[df_estuary['SAMPLE_No'].isin(needed_samples)][['SAMPLE_No', 'Easting', 'Northing']]

df_final = pd.concat([df_final,df_filtered])
print(len(df_final['SAMPLE_No']))