import numpy as np 
import matplotlib.pyplot as plt 
import sys 
import pyana 
import matplotlib

input_fitted = sys.argv[1] #fitted cube
input_obs = sys.argv[2]
input_lambda = sys.argv[3] #lambda grid
input_nodes = sys.argv[4]

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

for i in range(0,0):
	for j in range (0,0):
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
plt.cla()		

#now let's plot some "maps"

parameters = np.loadtxt(input_nodes,unpack=True)

temp = parameters.shape
NN = temp[0]

NX = 15 
NY = 15 #hardcoded

parameters = parameters.reshape(NN,NX,NY)

parameters[7] /= 1E5 #convert to km/s

matplotlib.rcParams['figure.figsize'] = 16,9
plt.subplot(321)
plt.pcolormesh(parameters[2])
plt.colorbar()
plt.title("Temperature at $\\log\\tau=-5$")
plt.subplot(322)
plt.pcolormesh(parameters[3])
plt.colorbar()
plt.title("Temperature at $\\log\\tau=-3.2$")
plt.subplot(323)
plt.pcolormesh(parameters[4])
plt.colorbar()
plt.title("Temperature at $\\log\\tau=-1.4$")
plt.subplot(324)
plt.pcolormesh(parameters[5])
plt.colorbar()
plt.title("Temperature at $\\log\\tau=0.5$")
plt.subplot(325)
plt.pcolormesh(parameters[7])
plt.colorbar()
plt.title("Velocity$")
plt.subplot(326)
plt.pcolormesh(obs_cube[0:15,0:15,0,30])
#plt.pcolormesh(parameters[8])
plt.colorbar()
plt.title("Continuum intensity")
plt.tight_layout()
plt.savefig("maps_test.png")
plt.clf()
plt.cla()

#plt.subplot(211)

plt.savefig("map_observed.png")








