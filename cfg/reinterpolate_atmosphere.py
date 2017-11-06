import numpy as np 
import matplotlib.pyplot as plt 
import scipy.interpolate as interpolate
import sys

input_atmosphere = sys.argv[1]
output_atmosphere = sys.argv[2]
ND = int(sys.argv[3])
index = int(sys.argv[4]) # which index we use to interpolate the atmosphere 
                         # typically it is 0 or 1 (tau or height) 

atmos_in = np.loadtxt(input_atmosphere, skiprows = 1)

dims = atmos_in.shape
ND_old = dims[0]
N_param = dims[1]

print ND_old, N_param

atmos_out = np.zeros((ND, N_param))

#start by making an independent variable, which is, of course, h
if (index == 0):
	atmos_out[:,0] = np.linspace(atmos_in[0,0], atmos_in[-1,0], num = ND)
	atmos_out[:,0] = np.linspace(-5.0, 1.2, num = ND)
elif (index == 1):
	atmos_out[:,1] = np.linspace(-atmos_in[0,1], -atmos_in[-1,1], num = ND)

for i in range (0,N_param):
	if (index == 0):
		f = interpolate.interp1d(atmos_in[:,0], atmos_in[:,i])
		atmos_out[:,i] = f(atmos_out[:,0])

	elif (index == 1):
		f = interpolate.interp1d(-atmos_in[:,1], atmos_in[:,i])
		atmos_out[:,i] = f(atmos_out[:,1])

if (index == 1):
	atmos_out[:,1] *= -1.0

np.savetxt(output_atmosphere, atmos_out, fmt = "%6.6e", header = str(ND) + " " + output_atmosphere)	