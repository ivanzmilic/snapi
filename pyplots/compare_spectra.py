import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage.filters as flt
import sys

file1 = sys.argv[1]
file2 = sys.argv[2]
output = sys.argv[3]
v_macro = float(sys.argv[4])

spectra1 = np.loadtxt(file1, unpack = True)
spectra2 = np.loadtxt(file2, unpack = True)

spectra1[1] /= max(spectra1[1])
spectra2[1] /= max(spectra2[1])

if spectra1[0][0] < 1:
	spectra1[0] *= 1E8
if spectra2[0][0] < 1:
	spectra2[0] *= 1E8

#spectra2[1,:] *= 2.997E29 / spectra2[0,:] / spectra2[0,:]

#print spectra2[1,:]

sigma = v_macro #km/s
sigma *= (spectra1[0,-1] + spectra1[0,0]) * 0.5 / 2.997E5 # Angstrom
sigma /= (spectra1[0,1]-spectra1[0,0])


spectra1[1] = flt.gaussian_filter(spectra1[1],sigma)

lambda_min = min([spectra1[0,0],spectra2[0,0]])
lambda_max = max([spectra1[0,-1],spectra2[0,-1]])

lambda_min = spectra1[0,0]
lambda_max = spectra1[0,-1]

plt.plot(spectra1[0,:], spectra1[1,:], color = 'red', label = file1)
plt.plot(spectra2[0,:], spectra2[1,:],color = 'blue', label = file2)
plt.xlim([lambda_min,lambda_max])
plt.ylim([0,1.2])
#plt.ylim([0,1E14])
plt.xlabel("$\lambda [\AA]$")
plt.ylabel("Normalized Intensity")
plt.legend()
plt.savefig(output)