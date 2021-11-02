import pandas as pd
import numpy as np

df = pd.read_csv('DATA/converted_chem_data.csv')
df = np.log10(df)
variance = df.var()
print(variance.sort_values(ascending=False))