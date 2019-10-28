from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)

import matplotlib
matplotlib.use('Agg')
import pyana
import numpy as np 
import matplotlib.pyplot as plt 
import sys
from scipy.signal import argrelextrema
import scipy.ndimage.filters as flt
from matplotlib import ticker 
from astropy.io import fits


cube1_in = sys.argv[1]
cube2_in = sys.argv[2]
stn = float(sys.argv[3])
ifmask = sys.argv[4]
maskfile = sys.argv[5]
fmt = sys.argv[6]

if (fmt == 1):
	temp = pyana.fzread(cube1_in)
	cube1 = temp["data"]
	temp = pyana.fzread(cube2_in)
	cube2 = temp["data"]
else:
	cube1 = fits.open(cube1_in)[0].data
	cube2 = fits.open(cube2_in)[0].data

NL = cube1.shape[-1]

if (int(ifmask)):
	mask = np.loadtxt(maskfile,skiprows=1)
	print mask.shape
	mask = mask[:NL]


cube_1_mean = np.mean(cube1[:,:,0,:],axis=(0,1))

noise = stn * np.sqrt(cube_1_mean[0]*cube_1_mean[:])
cube1 -= cube2
residual = cube1
del cube2
cube2 = 0.0 #wierd way to clear memory! 

weights = [float(sys.argv[6]),float(sys.argv[7]),float(sys.argv[8]),float(sys.argv[9])]
for i in range(0,4):
	weights[i] *= weights[i]
residual[:] /= noise;
residual *= residual
print residual.shape
if (int(ifmask)):
	residual[:,:,:,:] *= mask
residual = np.sum(residual,axis=3)
residual *= weights
residual = np.sum(residual,axis=2)

plt.clf()
plt.cla()
plt.imshow(np.log(residual))
plt.colorbar()
plt.savefig('chisq_map.png')

print 'chisq_max = ', np.amax(residual)
print 'chisq_mean = ', np.mean(residual)
print 'chisq_median = ', np.median(residual)
print weights
weights = np.asarray(weights)
N_Stokes = float(len(np.asarray(np.where(weights > 0.0)[0])))
if (int(ifmask)):
	residual /= N_Stokes * np.sum(mask)
else:
	residual /= (N_Stokes*NL)
print 'chisq_reduced_max = ', np.amax(residual)
print 'chisq_reduced_mean = ', np.mean(residual)
print 'chisq_reduced_median = ', np.median(residual)
