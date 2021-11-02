import numpy as np
import matplotlib.pyplot as plt

elem = 'Sr' #change to correct element
input_path = 'DATA/INVERSE_RESULTS/' + elem + '_results/' + elem + '_all_roughness_misfit.txt'
output_path = 'DATA/INVERSE_RESULTS/' + elem + '_tradeoff_plot.png'
input_array = np.loadtxt(input_path)

#setting up figure:
fig = plt.figure(figsize=(16,12))

#plotting misfit by roughness:
plt.scatter(input_array[3], input_array[1]+ input_array[2], c=input_array[0], cmap = 'viridis')

#add annotations of lambda:
for i,l in enumerate(input_array[0]):
    plt.annotate(l,(input_array[3][i], input_array[1][i]+ input_array[2][i]))

#adding labels and colourbar:
cb = plt.colorbar()
cb.set_label('lambda in log10')
plt.xlabel('Misfit')
plt.ylabel('Roughness (x+y)')
plt.title(elem+ ' tradeoff curve')
plt.savefig(output_path)