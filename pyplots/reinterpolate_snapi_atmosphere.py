import numpy as np 
import matplotlib.pyplot as plt 
import scipy.interpolate as interpolate
from astropy.io import fits
import pyana
import sys

def reinterpolate_atm_cube(atmos_in, grid, gridtype=1): # so far only tau

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
	if (gridtype == 1):
		atmos_out[1,None, None,None,:] = grid[:]
	else:
		atmos_out[0,None, None,None,:] = grid[:]
		#exit()

	if (gridtype == 1): # Flip because z goes in reverse
		atmos_in[:,:,:,:] = atmos_in[:,:,:,::-1]
		
	# take log of pressure
	atmos_in[3,:,:,:] = np.log10(atmos_in[3,:,:,:])
	atmos_in[4,:,:,:] = np.log10(atmos_in[4,:,:,:])

	for p in  range(0,NP):
		print ("interpolating the parameter p = ", p)
		for i in range(0,NX):
			for j in range(0,NY):
				
				f = interpolate.interp1d(atmos_in[gridtype,i,j,:], atmos_in[p,i,j,:], fill_value='extrapolate', kind='cubic')
				atmos_out[p,i,j,:] = f(grid)

	atmos_out[3,:,:,:] = 10.**atmos_out[3,:,:,:]
	atmos_out[4,:,:,:] = 10.**atmos_out[4,:,:,:]

	if (gridtype == 1): # Flip because z goes in reverse
		atmos_out[:,:,:,:] = atmos_out[:,:,:,::-1]

	return atmos_out
	
	

input_atmosphere_file = sys.argv[1]
output_atmosphere_file = sys.argv[2]

if (input_atmosphere_file[-3:] == '.f0'):
	atmos_og = pyana.fzread(input_atmosphere_file)["data"]
elif (input_atmosphere_file[-3:] == 'its'):
	atmos_og = fits.open(input_atmosphere_file)[0].data
elif (input_atmosphere_file[-4:] == '.dat'):
	atmos_og = np.loadtxt(input_atmosphere_file,skiprows=1, unpack=True)
	# reshape the atmosphere to the right dimensions
	atmos_og = atmos_og.reshape([atmos_og.shape[0], 1,1, atmos_og.shape[-1]])
else:
	print("info::not gonna work - I don't know the format")
	exit();

print("info::read the atmosphere with the dimensions: ", atmos_og.shape)


NDnew = int(sys.argv[3])
index = int(sys.argv[4]) # which index we use to interpolate the atmosphere 
                         # typically it is 0 or 1 (tau or height) 

lower = float(sys.argv[5])
upper = float(sys.argv[6])

grid_new = np.linspace(lower, upper, NDnew)

atmos_out = reinterpolate_atm_cube(atmos_og, grid_new, gridtype=index)

#kek = fits.PrimaryHDU(atmos_out)
#kek.writeto(output_atmosphere, overwrite=True)

# save as .dat just for the lolz
np.savetxt(output_atmosphere_file, atmos_out[:,0,0,:].T, header=str(atmos_out.shape[-1])+' OUT', fmt='%1.6e', comments='')
