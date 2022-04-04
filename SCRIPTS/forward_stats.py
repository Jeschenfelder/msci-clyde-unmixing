import numpy as np

elem = 'Zn' #<<<<<<<<<<<<<<<<<< change to correct element
lam = -0.2
#forward I/O:
#input_path = 'DATA/FORWARDMODEL_RESULTS/' + elem  + '_obs_v_pred.txt'
#output_path = 'DATA/FORWARDMODEL_RESULTS/' + elem + '_R2_misfit.txt'
#inverse I/O:
input_path = 'DATA/INVERSE_RESULTS/' + elem + '_results/' + elem +'_obs_v_pred.txt'
output_path = 'DATA/INVERSE_RESULTS/' + elem + '_results/' + elem + '_R2_misfit.txt'
#load in model observations vs predictions:
input_array = np.loadtxt(input_path,skiprows=1).T #6D array giving Sample_no X Y C_obs C_pred Misfit


#calculate rms misfit:
rms = np.sqrt(np.nansum(input_array[5]**2)/len(input_array[5]))

#calculate R2:
filtered_obs = input_array[4][~np.isnan(input_array[4])] #filtering observations to remove nan
filtered_pred = input_array[3][~np.isnan(input_array[4])] #filtering predictions to remove points where obs is nan
correlation_matrix = np.corrcoef(filtered_obs, filtered_pred)
R2 = correlation_matrix[0,1]**2

#filter out bad estuary points:
mask = np.logical_or(input_array[0]<700000, input_array[0]>700005)
masked_array = input_array[:,mask]

#calculate rms misfit:
rms_masked = np.sqrt(np.nansum(masked_array[5]**2)/len(masked_array[5]))

#calculate R2:
filtered_obs = masked_array[4][~np.isnan(masked_array[4])] #filtering observations to remove nan
filtered_pred = masked_array[3][~np.isnan(masked_array[4])] #filtering predictions to remove points where obs is nan
correlation_matrix = np.corrcoef(filtered_obs, filtered_pred)
R2_masked = correlation_matrix[0,1]**2

#save results:
out_array = np.array([rms,R2,rms_masked,R2_masked])


np.savetxt(output_path, out_array, fmt='%.3f')
