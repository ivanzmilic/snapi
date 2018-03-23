import pyana
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.ndimage.filters as flt
import scipy.interpolate as interpol
import sys

file_in = sys.argv[1]
file_out = sys.argv[2]

to_degrade = int(sys.argv[3])

temp = pyana.fzread(file_in)
stokes_cube = temp["data"]

dims = stokes_cube.shape
NX = dims[0]
NY = dims[1]
NL = 1501


#hardcode wavelength, in principle we could also read it (just a thought)
# in Angstrom
#l = np.linspace(8540.0,8543.0,NL)
l = np.linspace(15640.0,15670.0,NL)

sigma = -5 #in mA
sigma /= 2.35
sigma *= 1E-3  / (l[1]-l[0]) #to convert in 'pixels'
print sigma
noise_level = 1E-5
I_c_mean = np.mean(stokes_cube[:,:,0,0])
noise_level *= I_c_mean

plt.clf()
plt.cla()
plt.plot(l,stokes_cube[0,0,0],color='red')

#if we want to smear spatially:
if (to_degrade):
	A=[0.7,0.3] #two part-psf, weights
	width = [0.25,5.0] # two part - psf, widths in "
	A = np.asarray(A)
	width = np.asarray(width)
	width *= 725.0
	width /= 20.8
	width /= 2.35 #from FWHM to sigma
	print width

	for s in range (0,4):
		for w in range(0,NL):
			stokes_cube[:,:,s,w] = A[0] * flt.gaussian_filter(stokes_cube[:,:,s,w],width[0]) + A[1] * flt.gaussian_filter(stokes_cube[:,:,s,w],width[1])

print stokes_cube.shape

for i in range(0,NX):
	for j in range(0,NY):
		loc_noise = noise_level * np.sqrt(stokes_cube[i,j,0,0]/I_c_mean)
		for s in range(0,4):
			if (sigma>0):
				stokes_cube[i,j,s] = flt.gaussian_filter(stokes_cube[i,j,s],sigma)
			random_sample = np.random.normal(0,1.0,NL)
			stokes_cube[i,j,s] += random_sample*loc_noise

#Then we need to resample
NL_new = 601
#NL_new = 151
l_new = np.linspace(15640.0,15670.0,NL_new)
#l_new = np.linspace(8540.0,8543.0,NL_new)

print 'New sampling = ', l_new[1]-l_new[0]

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



