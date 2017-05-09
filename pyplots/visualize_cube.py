import numpy as np 
import matplotlib.pyplot as plt 
import sys 
import pyana 
import matplotlib

input_fitted = sys.argv[1] #fitted cube
input_lambdaf = sys.argv[2] #lambda used for fitting

input_obs = sys.argv[3] #observations we have fitted, not necessary same dimension as fitted
input_lambdao = sys.argv[4] #lambda of original observations, same as above

input_nodes = sys.argv[5] #inferred values of the nodes

matplotlib.rcParams['figure.figsize'] = 7, 10

#offsets between figures, hardcoded at the moment, we need to figure out how to put it in
x_offset = 5
y_offset = 5
l_offset = 677


#l = np.loadtxt(input_lambda)
#l*=1E8

a = pyana.fzread(input_fitted)
fitted_cube = a["data"]
dims = fitted_cube.shape
nx = dims[0]
ny = dims[1]
nlambda = dims[3]

b = pyana.fzread(input_obs)
obs_cube = b["data"]

print fitted_cube.shape
print obs_cube.shape

#These are debug lines
for i in range(0,5):
	for j in range (0,5):
		plt.subplot(211)
		plt.plot(fitted_cube[i,j,0,:290])
		plt.plot(obs_cube[j+x_offset,i+y_offset,0,677:])
		plt.savefig("test_cube"+str(i)+"_"+str(j)+".png")
		plt.xlabel("Wavelength")
		plt.ylabel("Stokes I")
		plt.subplot(212)
		#print fitted_cube[i,j,3,:290]
		plt.plot(fitted_cube[i,j,3,:290])
		plt.plot(obs_cube[j+x_offset,i+y_offset,3,677:])
		plt.savefig("test_cube"+str(i)+"_"+str(j)+".png")
		plt.xlabel("Wavelength")
		plt.ylabel("Stokes V")
		plt.tight_layout()
		plt.clf()
plt.cla()		

#now let's start by plotting observed map

#obsmap = plt.figure()
#obsmap.ion()
#obsmap.clf()
#obsmap.pcolormesh(obs_cube[:,:,0,0]) #plots stokes I at the first wavelength by default

#current_spectrum = plt.figure()
#current_spectrum.ion()
#current_spectrum.clf()
#current_spectrum.subplot(221)
#current_spectrum.plt(obs_cube[0,0,0,l_offset:l_offset+nlambda],'o')
#current_spectrum.subplot(222)
#current_spectrum.plt(obs_cube[0,0,1,l_offset:l_offset+nlambda],'o')
#current_spectrum.subplot(223)
#current_spectrum.plt(obs_cube[0,0,2,l_offset:l_offset+nlambda],'o')
#current_spectrum.subplot(224)
#current_spectrum.plt(obs_cube[0,0,3,l_offset:l_offset+nlambda],'o')
#current_spectrum.tight_layout()


parameters = np.loadtxt(input_nodes,unpack=True)

temp = parameters.shape
NN = temp[0]

NX = nx
NY = ny

parameters = parameters.reshape(NN,NX,NY)

parameters[7] /= 1E5 #convert to km/s

matplotlib.rcParams['figure.figsize'] = 16,9
plt.subplot(321)
plt.pcolormesh(parameters[2])
plt.ylim([0,NX])
plt.xlim([0,NY])
plt.colorbar()
plt.title("Temperature at $\\log\\tau=-5$")
plt.subplot(322)
plt.pcolormesh(parameters[3])
plt.ylim([0,NX])
plt.xlim([0,NY])
plt.colorbar()
plt.title("Temperature at $\\log\\tau=-3.2$")
plt.subplot(323)
plt.pcolormesh(parameters[4])
plt.ylim([0,NX])
plt.xlim([0,NY])
plt.colorbar()
plt.title("Temperature at $\\log\\tau=-1.4$")
plt.subplot(324)
plt.pcolormesh(parameters[5])
plt.ylim([0,NX])
plt.xlim([0,NY])
plt.colorbar()
plt.title("Temperature at $\\log\\tau=0.5$")
plt.subplot(325)
#plt.pcolormesh(parameters[7])
plt.pcolormesh(obs_cube[x_offset:x_offset+NX,y_offset:y_offset+NY,0,160+l_offset])
plt.ylim([0,NX])
plt.xlim([0,NY])
plt.colorbar()
plt.title("Line core intensity")
plt.subplot(326)
plt.pcolormesh(obs_cube[x_offset:x_offset+NX,y_offset:y_offset+NY,0,30])
plt.ylim([0,NX])
plt.xlim([0,NY])
#plt.pcolormesh(parameters[8])
plt.colorbar()
plt.title("Continuum intensity")
plt.tight_layout()
plt.savefig("maps_test.png")
plt.clf()
plt.cla()


#plt.subplot(211)

plt.savefig("map_observed.png")








