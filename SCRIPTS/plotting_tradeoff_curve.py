import numpy as np
import matplotlib.pyplot as plt

elem = 'Mg' #change to correct element
input_path = 'DATA/INVERSE_RESULTS/' + elem + '_results/' + elem + '_all_roughness_misfit.txt'
output_path = 'DATA/INVERSE_RESULTS/' + elem + '_tradeoff_plot.png'

#loading input and sorting by lambda:
input_array = np.loadtxt(input_path)

sorted_input = input_array[:, input_array[0].argsort()]
# setting misfit and calculating total roughness:
misfit = sorted_input[3]
roughness = sorted_input[1] + sorted_input[2] # total roughness as sum of x and y

#setting up figure:
fig = plt.figure(figsize=(16,12))

#plotting misfit by roughness:
plt.scatter(misfit, roughness
    , c=sorted_input[0], cmap = 'viridis')

#add annotations of lambda:
for i,l in enumerate(sorted_input[0]):
    plt.annotate(l,(misfit[i]+0.01, roughness[i]+0.005))

#adding labels and colourbar:
cb = plt.colorbar()
cb.set_label('lambda in log10')
plt.xlabel('Misfit')
plt.ylabel('Roughness x+y')
plt.title(elem + ' tradeoff curve')
plt.savefig(output_path)
#plt.show()