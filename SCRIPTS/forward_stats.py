import numpy as np

elem = 'Sr' #<<<<<<<<<<<<<<<<<< change to correct element
input_path = 'DATA/FORWARDMODEL_RESULTS/' + elem  + '_obs_v_pred.txt'
output_path = 'DATA/FORWARDMODEL_RESULTS/' + elem + '_R2_misfit.txt'
#load in model observations vs predictions:
input_array = np.loadtxt(input_path,skiprows=1).T #6D array giving Sample_no X Y C_obs C_pred Misfit

#calculate rms misfit:
rms = np.sqrt(np.nansum(input_array[5]**2)/len(input_array[5]))

#calculate R2:
filtered_obs = input_array[4][~np.isnan(input_array[4])] #filtering observations to remove nan
filtered_pred = input_array[3][~np.isnan(input_array[4])] #filtering predictions to remove points where obs is nan
correlation_matrix = np.corrcoef(filtered_obs, filtered_pred)
R2 = correlation_matrix[0,1]**2
#save results:
out_array = np.array([rms,R2])

np.savetxt(output_path, out_array, fmt='%.3f')
