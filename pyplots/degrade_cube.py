import pyana
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.ndimage.filters as flt
import scipy.interpolate as interpol
import sys

file_in = sys.argv[1]
file_out = sys.argv[2]

temp = pyana.fzread(file_in)
stokes_cube = temp["data"]

dims = stokes_cube.shape
NX = dims[0]
NY = dims[1]
NL = 1001


#hardcode wavelength, in principle we could also read it (just a thought)
# in Angstrom
l = np.linspace(5887.0,5897.0,NL)

sigma = 50 #in mA
sigma *= 1E-3  / (l[1]-l[0]) #to convert in 'pixels'
print sigma
noise_level = 3E-5
noise_level *= np.mean(stokes_cube[:,:,0,0])

plt.clf()
plt.cla()
plt.plot(l,stokes_cube[0,0,0],color='red')


if (sigma > 0):
	for i in range(0,NX):
		for j in range(0,NY):
			loc_noise = noise_level * np.sqrt(stokes_cube[i,j,0,:]/stokes_cube[i,j,0,20])
			for s in range(0,4):
				stokes_cube[i,j,s] = flt.gaussian_filter(stokes_cube[i,j,s],sigma)
				random_sample = np.random.normal(0,1.0,NL)
				stokes_cube[i,j,s] += random_sample*loc_noise

#Then we need to resample
NL_new = 501
l_new = np.linspace(5887.0,5897.0,NL_new)

resampled_cube = np.zeros([NX,NY,4,NL_new])

for i in range(0,NX):
	for j in range(0,NY):
		for s in range(0,4):
			f = interpol.interp1d(l,stokes_cube[i,j,s])
			resampled_cube[i,j,s] = f(l_new)

#After everything the Q,U,V need to be normalized:
for s in range(1,4):
	resampled_cube[:,:,s,:] /= resampled_cube[:,:,0,:]


plt.plot(l_new,resampled_cube[0,0,0],color='blue')
plt.savefig('pre_vs_post_degraded',fmt='png')

pyana.fzwrite(file_out,resampled_cube,0,'placeholder')



