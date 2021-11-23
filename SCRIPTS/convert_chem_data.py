#from os import error
#from numpy.lib.arraysetops import isin
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#defining conversion function:
def convert_to_format(df):
    cols_to_convert = []
    for col, data in df.iteritems(): # adding columns that need conversion from % to ppm to list for later use
        if data[0] == '%':
            cols_to_convert.append(col)
    df.drop(0, inplace=True) #drop first row
    for col in df.columns: #changing strings to NAN in columns
        df[col] = pd.to_numeric(df[col], errors='coerce')
        if col in cols_to_convert: #convert all concentrations to ppm (mg/kg)
            df[col] = df[col].multiply(10**4)
    return df

#defining oxide to elemental function:
def ox_to_elem(df):
    conversion_factors=np.array([['Na2O', 1.3480], ['MgO', 1.6582], ['Al2O3', 1.8895], ['SiO2', 2.1392], 
    ['P2O5', 2.2916], ['CaO', 1.3992], ['TiO2', 1.6681], ['Fe2O3', 1.4297], ['K2O', 1.2046]])
    for i in np.arange(0,len(conversion_factors)):
        oxide = conversion_factors[i][0]
        factor = conversion_factors[i][1].astype(float)
        df[oxide] = df[oxide].div(factor, fill_value=None)
    return(df)

################### work on CUSP Clyde Higher Order dataset: ###############################################
df_CUSP = pd.read_csv('DATA/ClydeHighOrderDrainage.dat',delimiter='\t')
#remove unnecessary rows (Easting, Northing, SAMPLE_TYPE, METHOD)
df_CUSP = df_CUSP.drop(['Easting', 'Northing', 'SAMPLE_TYPE', 'METHOD'], axis=1)
#remove rows with too much bad data (SO3, Cl, S, Cl, Ta, Tl, Ag, In, Te, Cd)
df_CUSP = df_CUSP.drop(['SO3', 'Cl', 'S', 'Cl.1', 'Ta', 'Tl', 'Ag', 'In', 'Te', 'Cd'], axis=1)
#print(df_CUSP.columns)

#make any bad data to NAN in cells:
df_CUSP = convert_to_format(df_CUSP)

# convert all oxides to elemental compositions:
df_CUSP = ox_to_elem(df_CUSP)


######################## working on Surface data set: #################################################################
df_surface = pd.read_csv('DATA/Compiled_Clyde_Surface.dat')
#remove unnecessary data: (SAMPLE_NAME, CRUISE, Lab_Number, Sample, Site, Easting, Northing, LOI, TOC, Hg)
df_surface = df_surface.drop(['SAMPLE_NAME', 'CRUISE', 'Lab_Number', 'Sample', 'Site', 'Easting', 'Northing', 'LOI', 'TOC'], axis=1)
#remove elements with too much bad data: (S, Cl, Ge, Se, Mo, Ag, Cd, In, Sb, Te, I, Cs, Yb, Ta, Ti, Bi, Hg)
df_surface = df_surface.drop(['S', 'Cl', 'Ge', 'Se', 'Mo', 'Ag', 'Cd', 'In', 'Sb', 'Te', 'I', 'Cs', 'Yb', 'Ta', 'Tl', 'Bi', 'Hg'], axis=1)

#make any bad data to NAN in cells:
df_surface = convert_to_format(df_surface)

# convert all oxides to elemental compositions:
df_surface = ox_to_elem(df_surface)

#filter out only needed samples:
sample_locs = np.loadtxt('DATA/filtered_sample_loc.dat')
sample_no = sample_locs.transpose()[2] #extract sample numbers from filtered data

for row, data in df_surface.iterrows():
    if data['SAMPLE_No'] in sample_no:
        continue
    else:
        df_surface.drop(row, inplace=True)
#adding sanity check that dataset is complete:
print('Samples in Surface data: ', len(df_surface.index))
print('Samples in CUSP data: ', len(df_CUSP.index))
############################################ concatinate and save: ###################################################################################
df_final = df_CUSP.copy() #create a 'final' dataframe

df_final = df_final.append(df_surface) #add surface data
print('Samples in final data set: ', len(df_final.index))
print(df_final.columns)

df_final.to_csv('DATA/converted_chem_data.csv', index=False)

'''
#plot duplicate columns: CaO, TiO2, Fe2O3, Ba
first = np.log10(df_CUSP['CaO'])
second = np.log10(df_CUSP['CaO.1'])
x = np.linspace(3.5,4.4, 100)
y = x
plt.plot(x,y,'k-')
plt.plot(first,second, 'rx')
plt.ylabel('log(mg/kg)')
plt.title('Duplicate CaO concentrations')
plt.show()
'''