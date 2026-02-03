import numpy as np 
import matplotlib.pyplot as plt
import sys

filein = sys.argv[1]
fileout = sys.argv[2]

atmos = np.loadtxt(filein, skiprows=1, unpack=True)

print("info:: shape of the input atmosphere:", atmos.shape)

# Here we will hard-code the magnetic field change.

B = 500.0
theta = 30.0 
phi = 0.0

atmos[7,:] = B
atmos[10,:] = np.radians(theta)
atmos[11,:] = np.radians(phi)

np.savetxt(fileout, atmos.T, fmt='%12.6e', header=str(atmos.shape[1])+' BLA', comments='')