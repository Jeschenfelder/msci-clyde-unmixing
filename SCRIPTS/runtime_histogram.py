import matplotlib.pyplot as plt
import numpy as np

times_array = np.loadtxt('msci-clyde-unmixing/DATA/timings.txt')/60**2
mean = np.mean(times_array)
plt.hist(times_array,color = 'green', label='Mean runtime (h): %.2f' %mean)
plt.ylabel('Count')
plt.xlabel('Runtime in hours')
plt.title('Histogram of runtimes for Mg inverse')
plt.yticks(np.arange(0,25,2))
plt.legend()
plt.show()