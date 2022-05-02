import numpy as np
import sys

def mass_per_km2(mass_low, mass_high,name,A):
    mass_avg = (mass_high+mass_low)/2
    m_km2_yr_low = mass_low/A
    m_km2_yr_high = mass_high/A
    m_km2_yr_avg = mass_avg/A
    print(name)
    print('Low: (kg km-2 yr-1) %.2e' %m_km2_yr_low)
    print('Average: (kg km-2 yr-1) %.2e' %m_km2_yr_avg)
    print('High: (kg km-2 yr-1 %.2e \n \n'%m_km2_yr_high)
    return

if __name__ == '__main__':
    r_earth = 6371 #km
    A_land = 4*np.pi*(r_earth**2)*0.3 #total land area

    name = sys.argv[1]
    m_low = float(sys.argv[2])*(10**6)
    m_high = float(sys.argv[3])*(10**6)

    mass_per_km2(m_low,m_high,name,A_land)
