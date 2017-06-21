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
input_atmos = sys.argv[6] #inferred atmospheres

print_maps_here = sys.argv[7]

#mask = np.loadtxt(sys.argv[7])

#matplotlib.rcParams['figure.figsize'] = 7, 10 #do I even want this?

#l_fit = np.loadtxt(input_lambdaf,unpack = True)
#l_obs = np.loadtxt(input_lambdao,unpack = True)

#l_fit[1]*=1E8
#l_obs[1]*=1E8

a = pyana.fzread(input_fitted)
fitted_cube = a["data"]
dims = fitted_cube.shape

l_offset = 495


#keep in mind this one is transposed:
NY = dims[0]
NX = dims[1]
NL = dims[3]

b = pyana.fzread(input_obs)
obs_cube = b["data"]

a_read = pyana.fzread(input_atmos)
atmospheres = a_read["data"]

print atmospheres.shape

#These are debug lines
for i in range(0,1):
	for j in range (0,1):
		plt.clf()
		plt.cla()
		plt.figure(figsize=[6,14])
		plt.subplot(311)
		plt.plot(fitted_cube[i,j,0,:])
		plt.plot(obs_cube[j,i,0,:])
		#plt.axvspan(669-l_offset,690-l_offset, alpha=0.5, color='red') #masks
		#plt.axvspan(725-l_offset,770-l_offset, alpha=0.5, color='red')
		#plt.axvspan(870-l_offset,905-l_offset, alpha=0.5, color='red')
		#plt.axvspan(915-l_offset,945-l_offset, alpha=0.5, color='red')
		plt.xlabel("Wavelength")
		plt.ylabel("Stokes I")
		plt.subplot(312)
		#print fitted_cube[i,j,3,:290]
		plt.plot(fitted_cube[i,j,3,:])
		plt.plot(obs_cube[j,i,3,:])
		plt.xlabel("Wavelength")
		plt.ylabel("Stokes V")

		plt.subplot(313)
		#print fitted_cube[i,j,3,:290]
		plt.plot(atmospheres[i,j,0],atmospheres[i,j,1])
		
		plt.xlabel("$\log \\tau$")
		plt.ylabel("$\mathrm{T\,[K]}$")
		plt.tight_layout()
		plt.savefig("test_cube"+str(i)+"_"+str(j)+".png")


#Here we read the parameter map.
pin = pyana.fzread(input_nodes)
parameters = pin["data"]
temp = parameters.shape
NN = temp[0] #total number of nodes
parameters[-2] /= 1E5 #convert to km/s										  
parameters[-1] /= 1E5 #convert to km/s


		
#print parameters[:,0,0]

if (NX<=1 and NY<=1):
	quit();




#Hard-coded plotting of some images.

barshrink = 0.8

plt.clf()
plt.cla()

#-----------------------------------------------------------------------------------------------------
# NOW NODES THEMSELVES -------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
panelsx=8
panelsy=2

#pre-determined wavelengths
l_core_Na = 342
l_core_Fe = 10
l_core_Ni = 28
l_c       = 140

plt.figure(figsize=[5*panelsx, 5*panelsy])

#Nodes panels:

plt.subplot(panelsy,panelsx,1)
plt.imshow(parameters[0],origin='lower')
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Temperature at $\log\,\\tau = -4$')

plt.subplot(panelsy,panelsx,2)
plt.imshow(parameters[1],origin='lower')
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Temperature at $\log \\tau = -3$')

plt.subplot(panelsy,panelsx,3)
plt.imshow(parameters[2],origin='lower')
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Temperature at $\log\,\\tau = -2$')

plt.subplot(panelsy,panelsx,4)
plt.imshow(parameters[3],origin='lower')
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Temperature at $\log \\tau = -1$')

plt.subplot(panelsy,panelsx,5)
plt.imshow(parameters[4],origin='lower')
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Temperature at $\log \\tau = 0$')

plt.subplot(panelsy,panelsx,6)
plt.imshow(parameters[5],origin='lower')
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('$\mathrm{v_{t}\,[km/s]}$')

plt.subplot(panelsy,panelsx,7)
plt.imshow(parameters[6],origin='lower')
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('$\mathrm{v_{los}\,[km/s]}$')

# Right hand side panels:

#Continuum intensity
i_cont = obs_cube[:,:,0,l_c].transpose()
i_c_mean = np.mean(i_cont)
i_cont /= i_c_mean
sigma = np.std(i_cont)

plt.subplot(panelsy,panelsx,9)
plt.imshow(i_cont,origin='lower',vmin = 1.0-3*sigma,vmax = 1.0+3*sigma)
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Observed continuum intensity')

i_cont = fitted_cube[:,:,0,l_c]
i_cont /= i_c_mean

plt.subplot(panelsy,panelsx,10)
plt.imshow(i_cont,origin='lower',vmin = 1.0-3*sigma,vmax = 1.0+3*sigma)
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Fitted continuum intensity')

#Fe intensity
i_cont = obs_cube[:,:,0,l_core_Fe].transpose()
i_c_mean = np.mean(i_cont)
i_cont /= i_c_mean
sigma = np.std(i_cont)

plt.subplot(panelsy,panelsx,11)
plt.imshow(i_cont,origin='lower',vmin = 1.0-3*sigma,vmax = 1.0+3*sigma)
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Observed Fe line core')

i_cont = fitted_cube[:,:,0,l_core_Fe]
i_cont /= i_c_mean

plt.subplot(panelsy,panelsx,12)
plt.imshow(i_cont,origin='lower',vmin = 1.0-3*sigma,vmax = 1.0+3*sigma)
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Fitted Fe line core')

#Ni intensity
i_cont = obs_cube[:,:,0,l_core_Ni].transpose()
i_c_mean = np.mean(i_cont)
i_cont /= i_c_mean
sigma = np.std(i_cont)

plt.subplot(panelsy,panelsx,13)
plt.imshow(i_cont,origin='lower',vmin = 1.0-3*sigma,vmax = 1.0+3*sigma)
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Observed Ni line core')

i_cont = fitted_cube[:,:,0,l_core_Ni]
i_cont /= i_c_mean

plt.subplot(panelsy,panelsx,14)
plt.imshow(i_cont,origin='lower',vmin = 1.0-3*sigma,vmax = 1.0+3*sigma)
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Fitted Ni line core')

#Na intensity
i_cont = obs_cube[:,:,0,l_core_Na].transpose()
i_c_mean = np.mean(i_cont)
i_cont /= i_c_mean
sigma = np.std(i_cont)

plt.subplot(panelsy,panelsx,15)
plt.imshow(i_cont,origin='lower',vmin = 1.0-3*sigma,vmax = 1.0+3*sigma)
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Observed Na D1 line core')

i_cont = fitted_cube[:,:,0,l_core_Na]
i_cont /= i_c_mean

plt.subplot(panelsy,panelsx,16)
plt.imshow(i_cont,origin='lower',vmin = 1.0-3*sigma,vmax = 1.0+3*sigma)
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Fitted Na D1 line core')



plt.tight_layout()
plt.savefig(print_maps_here+'.eps',fmt='eps')
plt.savefig(print_maps_here+'.png',fmt='png')
plt.clf()
plt.cla()

#------------------------------------------------------------------------------------------------------
#PLOTTING STRATIFICATION ------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------

panelsx=2
panelsy=4

l_core = 175


plt.figure(figsize=[10,12])

plt.subplot(panelsy*100+panelsx*10+1)
plt.imshow(atmospheres[:,:,1,1],origin='lower')
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Temperature at $\log\,\\tau = -4.9$')

plt.subplot(panelsy*100+panelsx*10+3)
plt.imshow(atmospheres[:,:,1,5],origin='lower')
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Temperature at $\log \\tau = -3.9$')

plt.subplot(panelsy*100+panelsx*10+5)
plt.imshow(atmospheres[:,:,1,8],origin='lower')
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Temperature at $\log \\tau = -3$')

plt.subplot(panelsy*100+panelsx*10+7)
plt.imshow(atmospheres[:,:,1,12],origin='lower')
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Temperature at $\log \\tau = -1.8$')

#Continuum intensity
i_cont = obs_cube[:,:,0,30].transpose()
i_c_mean = np.mean(i_cont)
i_cont /= i_c_mean
sigma = np.std(i_cont)

plt.subplot(panelsy*100+panelsx*10+8)
plt.imshow(i_cont,origin='lower',vmin = 1.0-3*sigma,vmax = 1.0+3*sigma)
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Continuum intensity')

#Core intensity
i_core = obs_cube[:,:,0,l_core].transpose()
i_core_mean = np.mean(i_core)
i_core /= i_core_mean
sigma = np.std(i_cont)

plt.subplot(panelsy*100+panelsx*10+6)
plt.imshow(i_core,origin='lower',vmin=1.0-3*sigma,vmax=1.0+3*sigma)
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Observed Na D1 line core intensity')

plt.subplot(panelsy*100+panelsx*10+2)
plt.imshow(atmospheres[:,:,1,15],origin='lower')
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Temperature at $\log \\tau = -0.9$')

plt.subplot(panelsy*100+panelsx*10+4)
plt.imshow(atmospheres[:,:,1,18],origin='lower')
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Temperature at $\log \\tau = 0$')


plt.tight_layout()
plt.savefig(print_maps_here+'_temp_strat.eps',fmt='eps')
plt.savefig(print_maps_here+'_temp_strat.png',fmt='png')
plt.clf()
plt.cla()














