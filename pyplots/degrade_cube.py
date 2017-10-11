import pyana
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.ndimage.filters as flt
import sys

file_in = sys.argv[1]
file_out = sys.argv[2]

temp = pyana.fzread(file_in)
stokes_cube = temp["data"]

dims = stokes_cube.shape
NX = dims[0]
NY = dims[1]


#hardcode wavelength, in principle we could also read it (just a thought)
# in Angstrom
l = np.linspace(5893.0,5897.0,401)

sigma = 20.0 #in mA
sigma *= 1E-3  / (l[1]-l[0]) #to convert in 'pixels'
noise_level = 5E-4
noise_level *= np.mean(stokes_cube[:,:,0,20])

plt.clf()
plt.cla()
plt.plot(l,stokes_cube[0,0,3],color='red')

for i in range(0,NX):
	for j in range(0,NY):
		loc_noise = noise_level * np.sqrt(stokes_cube[i,j,0,:]/stokes_cube[i,j,0,20])
		for s in range(0,4):
			stokes_cube[i,j,s] = flt.gaussian_filter(stokes_cube[i,j,s],sigma)
			random_sample = np.random.normal(0,1.0,401)
			stokes_cube[i,j,s] += random_sample*loc_noise

plt.plot(l,stokes_cube[0,0,3],'o',color='blue')
plt.savefig('pre_vs_post_degraded',fmt='png')

pyana.fzwrite(file_out,stokes_cube,0,'placeholder')



