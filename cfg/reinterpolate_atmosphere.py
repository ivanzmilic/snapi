import numpy as np 
import matplotlib.pyplot as plt 
import scipy.interpolate as interpolate
import sys

input_atmosphere = sys.argv[1]
output_atmosphere = sys.argv[2]
ND = int(sys.argv[3])

atmos_in = np.loadtxt(input_atmosphere, skiprows = 1)

atmos_out = np.zeros((ND, 11))

#start by making an independent variable, which is, of course, h

atmos_out[:,1] = np.linspace(-atmos_in[0,1], -atmos_in[-1,1], num = ND)

#print atmos_out[:,1]

for i in range (0,11):
	if i != 1:
		f = interpolate.interp1d(-atmos_in[:,1], atmos_in[:,i])
		atmos_out[:,i] = f(atmos_out[:,1])

atmos_out[:,1] *= -1.0

np.savetxt(output_atmosphere, atmos_out, fmt = "%6.6e", header = str(ND) + " " + output_atmosphere)	