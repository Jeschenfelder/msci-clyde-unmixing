import numpy as np
import matplotlib.pyplot as plt

elem = 'Mg' #change to correct element
input_path = 'DATA/INVERSE_RESULTS/' + elem + ''
input_array = np.loadtxt(input_path)

plt.scatter(input_array[1], input_array[2], c=input_array[0])
plt.show()