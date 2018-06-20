import pyana
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.ndimage.filters as flt
import scipy.interpolate as interpol
import sys

file_in = sys.argv[1]
file_out = sys.argv[2]
to_normalize = int(sys.argv[3])


temp = pyana.fzread(file_in)
stokes_cube = temp["data"]

dims = stokes_cube.shape
NX = dims[0]
NY = dims[1]
NL = dims[3]

print NX,NY,NL

noise_level = 3E-4
I_c_mean = np.mean(stokes_cube[:,:,0,0])
noise_level *= I_c_mean

for i in range(0,NX):
	for j in range(0,NY):
		loc_noise = noise_level * np.sqrt(stokes_cube[i,j,0,0]/I_c_mean)
		for s in range(0,4):
			random_sample = np.random.normal(0,1.0,NL)
			stokes_cube[i,j,s] += random_sample*loc_noise

if(to_normalize):
	for s in range(1,4):
		stokes_cube[:,:,s,:] /= stokes_cube[:,:,0,:]

pyana.fzwrite(file_out,stokes_cube,0,'placeholder')
