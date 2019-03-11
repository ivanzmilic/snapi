import pyana 
import numpy as np 
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import sys

file_in = sys.argv[1]
file_out = sys.argv[2]

temp = pyana.fzread(file_in)
atmos = temp["data"]

N_new = 41
tau_uniform = np.linspace(-5,1,N_new)

dims = atmos.shape
NX = dims[0]
NY = dims[1]

atmos_new = np.zeros([NX,NY,12,N_new])

for i in range(0,NX):
	for j in range(0,NY):
		atmos_new[i,j,0,:] = tau_uniform[:]
		for p in range (1,12):
			f = interp1d(atmos[i,j,0],atmos[i,j,p])
			atmos_new[i,j,p] = f(tau_uniform)

pyana.fzwrite(file_out,atmos_new,0,'placeholder')
