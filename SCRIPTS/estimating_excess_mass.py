import numpy as np

out_file = open('DATA/excess_mass.txt', 'w')

# loading in data:
inverse_profile = np.loadtxt('DATA/INVERSE_RESULTS/Pb_results/Pb_pred_profile.txt', delimiter=' ')
inverse_conc = inverse_profile[:,1]

forward_profile = np.loadtxt('DATA/FORWARDMODEL_RESULTS/Pb_pred_profile.txt', delimiter=' ')
forward_conc = forward_profile[:,1]

flux_profile = np.loadtxt('DATA/sedimentflux_profile.txt', delimiter=' ')
flux_mass = flux_profile[:,1]

# converting into total mass:
inverse_mass = (np.power(10, inverse_conc))*(10**(-6))*flux_mass
forward_mass = (np.power(10, forward_conc))*(10**(-6))*flux_mass

# calcluating mass difference:
diff_mass = inverse_mass - forward_mass
out_file.write("Pb excess mass (kg): %.2f \n" %(np.sum(diff_mass)))

########################################################################################################################################

# loading in data:
inverse_profile = np.loadtxt('DATA/INVERSE_RESULTS/Cu_results/Cu_pred_profile.txt', delimiter=' ')
inverse_conc = inverse_profile[:,1]

forward_profile = np.loadtxt('DATA/FORWARDMODEL_RESULTS/Cu_pred_profile.txt', delimiter=' ')
forward_conc = forward_profile[:,1]

# converting into total mass:
inverse_mass = (np.power(10, inverse_conc))*(10**(-6))*flux_mass
forward_mass = (np.power(10, forward_conc))*(10**(-6))*flux_mass

# calcluating mass difference:
diff_mass = inverse_mass - forward_mass
out_file.write("Cu excess mass (kg): %.2f \n" %(np.sum(diff_mass)))

########################################################################################################################################

# loading in data:
inverse_profile = np.loadtxt('DATA/INVERSE_RESULTS/Zn_results/Zn_pred_profile.txt', delimiter=' ')
inverse_conc = inverse_profile[:,1]

forward_profile = np.loadtxt('DATA/FORWARDMODEL_RESULTS/Zn_pred_profile.txt', delimiter=' ')
forward_conc = forward_profile[:,1]

# converting into total mass:
inverse_mass = (np.power(10, inverse_conc))*(10**(-6))*flux_mass
forward_mass = (np.power(10, forward_conc))*(10**(-6))*flux_mass

# calcluating mass difference:
diff_mass = inverse_mass - forward_mass
out_file.write("Zn excess mass (kg): %.2f \n" %(np.sum(diff_mass)))
out_file.close()