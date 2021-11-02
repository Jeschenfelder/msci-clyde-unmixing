from os import sep
import pandas as pd

df = pd.read_csv('DATA/GBASESEDIMENTS/prepped_gbase_auto.dat', delimiter=' ')
print(df.head())
df = df.fillna('NAN')
print(df.head())
df.to_csv('DATA/GBASESEDIMENTS/prepped_gbase_auto_NAN.dat', sep=' ')