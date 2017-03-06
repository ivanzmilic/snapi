import numpy as np 
import matplotlib.pyplot as plt 

matrix = np.loadtxt("../response_matrix.txt", unpack = True)
matrix = matrix.T

plt.pcolormesh(matrix, rasterized = True)
plt.colorbar()
plt.tight_layout()
plt.savefig("response_matrix.eps", fmt = 'eps')