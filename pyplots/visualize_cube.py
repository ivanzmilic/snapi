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

print_maps_here = sys.argv[6]

#mask = np.loadtxt(sys.argv[7])

#matplotlib.rcParams['figure.figsize'] = 7, 10 #do I even want this?

#l = np.loadtxt(input_lambda)
#l*=1E8

a = pyana.fzread(input_fitted)
fitted_cube = a["data"]
dims = fitted_cube.shape

#print fitted_cube.shape
#print dims
#keep in mind this one is transposed:
NY = dims[0]
NX = dims[1]
NL = dims[3]

b = pyana.fzread(input_obs)
obs_cube = b["data"]

#print obs_cube.shape

l_offset = 660

#print fitted_cube[0,0,3,:]

#These are debug lines
for i in range(0,1):
	for j in range (0,1):
		plt.clf()
		plt.cla()
		plt.figure(figsize=[6,10])
		plt.subplot(211)
		plt.plot(fitted_cube[i,j,0,:])
		plt.plot(obs_cube[j,i,0,:])
		plt.axvspan(669-l_offset,690-l_offset, alpha=0.5, color='red')
		plt.axvspan(725-l_offset,770-l_offset, alpha=0.5, color='red')
		plt.axvspan(870-l_offset,905-l_offset, alpha=0.5, color='red')
		plt.axvspan(915-l_offset,945-l_offset, alpha=0.5, color='red')
		plt.xlabel("Wavelength")
		plt.ylabel("Stokes I")
		plt.subplot(212)
		#print fitted_cube[i,j,3,:290]
		plt.plot(fitted_cube[i,j,3,:])
		plt.plot(obs_cube[j,i,3,:])
		plt.xlabel("Wavelength")
		plt.ylabel("Stokes V")
		plt.tight_layout()
		plt.savefig("test_cube"+str(i)+"_"+str(j)+".png")
		

if (NX<=1 and NY<=1):
	quit();


#Here we read the parameter map.
parameters = np.loadtxt(input_nodes,unpack=True)
temp = parameters.shape
NN = temp[0] #total number of nodes

parameters = parameters.reshape(NN,NX,NY) #we do this without checking but we should 
										  #probably check to see if dimensions match
parameters[7] /= 1E5 #convert to km/s

#Hard-coded plotting of some images.

barshrink = 0.8

plt.clf()
plt.cla()

panelsx=2
panelsy=4

l_core = 175


plt.figure(figsize=[10,12])

plt.subplot(panelsy*100+panelsx*10+1)
plt.imshow(parameters[2].transpose(),origin='lower')
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Temperature at $\log\,\\tau = -5$')

plt.subplot(panelsy*100+panelsx*10+3)
plt.imshow(parameters[3].transpose(),origin='lower',vmin=5500,vmax=7000)
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Temperature at $\log \\tau = -3.2$')

plt.subplot(panelsy*100+panelsx*10+5)
plt.imshow(parameters[4].transpose(),origin='lower')
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Temperature at $\log\,\\tau = -1.4$')

plt.subplot(panelsy*100+panelsx*10+7)
plt.imshow(parameters[5].transpose(),origin='lower')
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Temperature at $\log \\tau = 0.5$')

i_cont = obs_cube[:,:,0,30].transpose()

i_c_mean = np.mean(i_cont)
i_cont /= i_c_mean

sigma = np.std(i_cont)
#print i_c_mean, sigma

plt.subplot(panelsy*100+panelsx*10+8)
plt.imshow(i_cont,origin='lower',vmin = 1.0-3*sigma,vmax = 1.0+3*sigma)
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Continuum intensity')

i_core = obs_cube[:,:,0,l_core].transpose()
i_core_mean = np.mean(i_core)
i_core /= i_core_mean
sigma = np.std(i_cont)

plt.subplot(panelsy*100+panelsx*10+6)
plt.imshow(i_core,origin='lower',vmin=1.0-3*sigma,vmax=1.0+3*sigma)
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Observed Na D1 line core intensity')

#parameters[8] *= np.cos(parameters[9]*np.pi/180.0)

#plt.subplot(panelsy*100+panelsx*10+2)
#plt.imshow(parameters[8].transpose(),origin='lower',vmin=-1000,vmax=1000)
#plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
#plt.title('$\mathrm{B\,[Gauss]}$')

i_core = fitted_cube[:,:,0,l_core]
i_core_mean = np.mean(i_core)
i_core /= i_core_mean
sigma = np.std(i_cont)
plt.subplot(panelsy*100+panelsx*10+2)
plt.imshow(i_core,origin='lower',vmin=1.0-3*sigma,vmax=1.0+3*sigma)
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Fitted Na D1 line core intensity')



plt.subplot(panelsy*100+panelsx*10+4)
plt.imshow(parameters[7].transpose(),origin='lower')
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('$\mathrm{v_{los}\,[km/s]}$')

plt.tight_layout()
plt.savefig(print_maps_here+'.eps',fmt='eps')
plt.savefig(print_maps_here+'.png',fmt='png')
plt.clf()
plt.cla()











