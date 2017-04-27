import numpy as np 
import matplotlib.pyplot as plt 
import sys 
import pyana 
import matplotlib

input_fitted = sys.argv[1] #fitted cube
input_obs = sys.argv[2]
input_lambda = sys.argv[3] #lambda grid

matplotlib.rcParams['figure.figsize'] = 7, 10

l = np.loadtxt(input_lambda)
l*=1E8

a = pyana.fzread(input_fitted)
fitted_cube = a["data"]
dims = fitted_cube.shape
nx = dims[0]
ny = dims[1]

b = pyana.fzread(input_obs)
obs_cube = b["data"]

for i in range(0,nx):
	for j in range (14,15):
		plt.subplot(211)
		plt.plot(l,fitted_cube[i,j,0,:])
		plt.plot(l,obs_cube[i,j,0,677:])
		plt.savefig("test_cube"+str(i)+"_"+str(j)+".png")
		plt.xlabel("Wavelength")
		plt.ylabel("Stokes I")
		plt.subplot(212)
		plt.plot(l,fitted_cube[i,j,3,:])
		plt.plot(l,obs_cube[i,j,3,677:])
		plt.savefig("test_cube"+str(i)+"_"+str(j)+".png")
		plt.xlabel("Wavelength")
		plt.ylabel("Stokes V")
		plt.tight_layout()
		plt.clf()		


