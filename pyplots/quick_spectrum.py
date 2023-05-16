import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage.filters as flt
import sys

spectrum = np.loadtxt(sys.argv[1])
spectrum[:,0] *= 1E8

plt.clf()
plt.cla()
plt.figure(figsize=[6.5,3.5])


plt.plot(spectrum[:,0], spectrum[:,1])
lambda_l = min(spectrum[:,0])
lambda_m = max(spectrum[:,0])
#spectrum[:,1] /= max(spectrum[:,1])

#Other spectra
dims = spectrum.shape
n_spectra = dims[1] - 1;

#for i in range(2,n_spectra+1):
#	plt.plot(spectrum[:,0], spectrum[:,i])

#lambda_l = 2785.0
#lambda_m = 2810.0

#plt.ion()
plt.xlabel("$\mathrm{Wavelength [\AA]}$")
plt.ylabel("$\mathrm{Intensity}$")
plt.xlim([lambda_l, lambda_m])
#plt.ylim([0, np.max(spectrum[:,1])*1.1])
#plt.ylim([1,1.5E14])
plt.tight_layout()
plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
plt.savefig(sys.argv[2]+'.eps', format='eps',bbox_inches='tight')
plt.savefig(sys.argv[2]+'.png', bbox_inches='tight')

