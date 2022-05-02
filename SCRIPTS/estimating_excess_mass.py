import numpy as np

def excess_calc(elem):
    in_str = 'DATA/INVERSE_RESULTS/'+elem+'_results/'+elem+'_pred_profile.txt'
    for_str = 'DATA/FORWARDMODEL_RESULTS/'+elem+'_pred_profile.txt'

    # loading in data:
    inverse_profile = np.loadtxt(in_str, delimiter=' ')
    inverse_conc = inverse_profile[:,1]
    forward_profile = np.loadtxt(for_str, delimiter=' ')
    forward_conc = forward_profile[:,1]
    flux_profile = np.loadtxt('DATA/sedimentflux_profile.txt', delimiter=' ')
    flux_mass = flux_profile[:,1]

    # converting into total mass:
    inverse_mass = (np.power(10, inverse_conc))*(10**(-6))*flux_mass
    forward_mass = (np.power(10, forward_conc))*(10**(-6))*flux_mass
    print('Total Inverse mass of %s (kg/yr): ' %(elem),inverse_mass[0])
    print('Total natural mass of %s (kg/yr):' %elem, forward_mass[0])

    # calcluating mass difference:
    diff_mass = inverse_mass - forward_mass
    mass_profile = np.array([flux_profile[:,0],diff_mass]).T

    # calculating dM/dx:
    dM_dx = np.diff(diff_mass)/np.diff(flux_profile[:,0])
    dM_dx_profile = np.array([flux_profile[1:,0],dM_dx]).T


    #writing out total excess and as function of distance:
    out_str_1 = 'DATA/'+elem+'_excess_mass_profile.txt'
    out_str_2 = 'DATA/'+elem+'_mass_flux_profile.txt'
    np.savetxt(out_str_2,dM_dx_profile)
    np.savetxt(out_str_1,mass_profile)
    out_file.write("%s excess mass (kg): %.2f \n" %(elem, diff_mass[0]))

out_file = open('DATA/excess_mass.txt', 'w')

excess_calc('Pb')
excess_calc('Zn')
excess_calc('Cu')

out_file.close()