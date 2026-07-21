import numpy as np 
import matplotlib.pyplot as plt 
import scipy.interpolate as interpolate
from astropy.io import fits
import pyana
import sys

def reinterpolate_atm_cube(atmos_in, grid): # so far only tau

	dims = atmos_in.shape
	# implying orientation NP, NX, NY, NZ:
	NZ_old = dims[-1]
	NP = dims[0]
	NX = dims[1]
	NY = dims[2]

	print ("info::original number of depths and number of parameters are: ", NZ_old, NP)

	NZ = len(grid)

	atmos_out = np.zeros([NP, NX, NY, NZ])

	#start by making an independent variable, which is, gonna be tau
	atmos_out[None,0,None,None,:] = grid[:]


	# take log of pressure
	atmos_in[3,:,:,:] = np.log10(atmos_in[3,:,:,:])
	atmos_in[4,:,:,:] = np.log10(atmos_in[4,:,:,:])

	for p in  range(1,NP):
		print ("interpolating the parameter p = ", p)
		for i in range(0,NX):
			for j in range(0,NY):
				
				# it's a snapi atmosphere so logtau is increasing (top is at the start):
				f = interpolate.interp1d(atmos_in[0,i,j,:], atmos_in[p,i,j,:], fill_value='extrapolate', kind='cubic')
				atmos_out[p,i,j,:] = f(grid)

	atmos_out[3,:,:,:] = 10.**atmos_out[3,:,:,:]
	atmos_out[4,:,:,:] = 10.**atmos_out[4,:,:,:]
	return atmos_out
	
	

input_atmosphere = sys.argv[1]
output_atmosphere = sys.argv[2]

if (input_atmosphere[-3:] == '.f0'):
	atmos_og = pyana.fzread(input_atmosphere)["data"]
elif (input_atmosphere[-3:] == 'its'):
	atmos_og = fits.open(input_atmosphere)[0].data
else:
	print("info::not gonna work - I don't know the format")
	exit();

print("info::read the atmosphere with the dimensions: ", atmos_og.shape)


ND = int(sys.argv[3])
index = int(sys.argv[4]) # which index we use to interpolate the atmosphere 
                         # typically it is 0 or 1 (tau or height) 

tau_new = np.linspace(-4.0,1.0,ND)

atmos_out = reinterpolate_atm_cube(atmos_og, tau_new) 

kek = fits.PrimaryHDU(atmos_out)
kek.writeto(output_atmosphere, overwrite=True)
