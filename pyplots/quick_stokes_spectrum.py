import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage.filters as flt
import sys

spectrum = np.loadtxt(sys.argv[1])

scalar_spectrum = np.loadtxt(sys.argv[2])

if (spectrum[0,0] < 1.0):
	spectrum[:,0] *= 1E8
if (scalar_spectrum[0,0] < 1.0):
	scalar_spectrum[:,0] *=1E8
lambda_l = -1.5
lambda_m = 1.5
#spectrum[:,1] /= max(spectrum[:,1])

#To compare with fits atlas
#V_over_I = spectrum[:,4] / max(spectrum[:,1])
#V_over_I /= max(V_over_I/2.0)
#sigma = 2.0
#sigma *= 5896.0 / 2.997E5 # Angstrom
#sigma /= (spectrum[1,0]-spectrum[0,0])
#V_over_I = flt.gaussian_filter(V_over_I,sigma)

#plt.clf()
#plt.cla()
#plt.plot(spectrum[:,0],V_over_I)
#plt.plot(scalar_spectrum[:,0],scalar_spectrum[:,2])
#plt.xlim([5895,5897])
#plt.savefig("Na_comparison_with_fts_Stokes.png", format = 'png')


#We can compute the weak field approximation story:

V_weak_field = -np.gradient(spectrum[:,1]) / np.gradient(spectrum[:,0]) * 4.697E-13 * ((lambda_l + lambda_m) ** 2.0) * 0.25 * 2.0 * 1000.0

#spectrum[:,0] -= (spectrum[0,0] + spectrum[-1,0]) * 0.5

lambda_l = spectrum[0,0]
lambda_m = spectrum[-1,0]

plt.figure(1)
plt.clf()
plt.cla()
plt.subplot(221)
#plt.ion()
plt.plot(spectrum[:,0], spectrum[:,1])
if (sys.argv[1] != sys.argv[2]):
	plt.plot(spectrum[:,0], scalar_spectrum[:,1])
#plt.xlabel("Wavelength")
plt.ylabel("Intensity")
plt.xlim([lambda_l, lambda_m])
#plt.xlim([6301.4, 6301.6])
plt.ylim([0, np.max(spectrum[:,1])*1.1])
plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
plt.tight_layout()

plt.subplot(222)
#plt.ion()
plt.plot(spectrum[:,0], spectrum[:,2])
if (sys.argv[1] != sys.argv[2]):
	plt.plot(spectrum[:,0], scalar_spectrum[:,2])
#plt.xlabel("Wavelength")
#plt.ylabel("Intensity")
plt.xlim([lambda_l, lambda_m])
#plt.ylim([0, np.max(spectrum[:,1])*1.1])
plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
plt.tight_layout()

plt.subplot(223)
#plt.ion()
plt.plot(spectrum[:,0], spectrum[:,3])
if (sys.argv[1] != sys.argv[2]):
	plt.plot(spectrum[:,0], scalar_spectrum[:,3])
plt.xlabel("Wavelength")
plt.ylabel("Intensity")
plt.xlim([lambda_l, lambda_m])
#plt.ylim([0, np.max(spectrum[:,1])*1.1])
plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
plt.tight_layout()

plt.subplot(224)
#plt.ion()
plt.plot(spectrum[:,0], spectrum[:,4])
if (sys.argv[1] != sys.argv[2]):
	plt.plot(spectrum[:,0], scalar_spectrum[:,4])

#plt.plot(spectrum[:,0], V_weak_field)
plt.xlabel("Wavelength")
#plt.ylabel("Intensity")
plt.xlim([lambda_l, lambda_m])
#plt.ylim([, np.max(spectrum[:,4])*1.1])
plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
plt.tight_layout()

plt.savefig(sys.argv[3], format='png')
	
