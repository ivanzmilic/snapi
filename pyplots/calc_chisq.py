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

cube1_in = sys.argv[1]
cube2_in = sys.argv[2]
stn = float(sys.argv[3])
ifmask = sys.argv[4]
maskfile = sys.argv[5]

temp = pyana.fzread(cube1_in)
cube1 = temp["data"]
temp = pyana.fzread(cube2_in)
cube2 = temp["data"]

NL = cube1.shape[-1]

if (int(ifmask)):
	mask = np.loadtxt(maskfile,skiprows=1)
	print mask.shape

cube_1_mean = np.mean(cube1[:,:,0,:],axis=(0,1))

noise = stn * np.sqrt(cube_1_mean[0]*cube_1_mean[:])
cube1 -= cube2
residual = cube1
del cube2
cube2 = 0.0 #wierd way to clear memory! 

weigths = [1.0,0.0,0.0,2.0]
residual[:] /= noise;
residual *= residual
print residual.shape
if (int(ifmask)):
	residual[:,:,:,:] *= mask
residual = np.sum(residual,axis=3)
residual *= weigths**2.0
residual = np.sum(residual,axis=2)

plt.clf()
plt.cla()
plt.imshow(np.log(residual))
plt.colorbar()
plt.savefig('chisq_map.png')

print 'chisq_max = ', np.amax(residual)
print 'chisq_mean = ', np.mean(residual)
if (int(ifmask)):
	residual /= 2.0 * np.sum(mask)
else:
	residual /= (2.0*NL)
print 'chisq_reduced_max = ', np.amax(residual)
print 'chisq_reduced_mean = ', np.mean(residual)

