import matplotlib
matplotlib.use('Agg')
import numpy as np 
import scipy.signal as sgn
import matplotlib.pyplot as plt 
import sys 
import pyana 
import matplotlib
from scipy.signal import argrelextrema
import scipy.ndimage.filters as flt

# Written by I.Milic (MPS), 2017
# Very simple, and mostly hardcoded procedure to visualize the results of inversion. 
# Goal one: visualize the resulting atmospheres/nodes. 
# Goal two: compare fitted and observed data

#input_fitted = sys.argv[1] #fitted cube
#input_lambdaf = sys.argv[2] #lambda used for fitting

#input_obs = sys.argv[3] #observations we have fitted, not necessary same dimension as fitted
#input_lambdao = sys.argv[4] #lambda of original observations, same as above

input_nodes = sys.argv[1] #inferred values of the nodes
input_atmos = sys.argv[2] #inferred atmospheres
input_spectra = sys.argv[3]


print_maps_here = sys.argv[4]

#Here we read the parameter/node map:
temp = pyana.fzread(input_nodes)
parameters = temp["data"]
#parameters = np.transpose(parameters,(0,2,1))
temp = parameters.shape
NN = temp[0] #total number of nodes

a_read = pyana.fzread(input_atmos)
atmospheres = a_read["data"]

temp = pyana.fzread(input_spectra)
stokes = temp["data"]
stokes = np.transpose(stokes,(1,0,2,3))
parameters = np.transpose(parameters,(0,1,2))

I_c = np.mean(stokes[:,:,0,0])
stokes[:,:,:,:] /= I_c

#-----------------------------------------------------------------------------------------------------
# NOW PLOT THE NODES ---------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------


T_nodes_tau = [-3.0,-2.0,-1.0,0.0]
#T_nodes = np.arange(len(T_nodes_tau))
T_nodes = [1,3]
vs_nodes_tau = [-2.5,-1.0,0.0]
#vs_nodes = np.arange(len(vs_nodes_tau)) + len(T_nodes_tau)
vs_nodes = [5,7]
B_nodes_tau = [-2.0,-0.0]
#B_nodes = np.arange(len(B_nodes_tau)) + len(vs_nodes_tau) + len(T_nodes_tau)
B_nodes = [8,9]
theta_nodes = [len(B_nodes_tau) + len(vs_nodes_tau) + len(T_nodes_tau)]
theta_nodes=[10]

panelsx=len(T_nodes)
panelsy=4

Tmap = 'hot'
Vmap = 'coolwarm'
Bmap = 'cividis'
Imap = 'hot'
Pmap = 'Spectral'
Dmap = 'coolwarm'

plt.clf()
plt.cla()

panel_no = 1

defsize = 8
ratio = 0.5
barshrink=0.8
plt.figure(figsize=[defsize*panelsx, ratio*defsize*panelsy])
k =0

plt.subplot(panelsy,panelsx,panel_no)
plt.imshow(stokes[:,:,0,0],origin='lower',cmap=Imap,vmin=0.7,vmax=1.3)
#if (i==2):
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Stokes $I$ (continuum)')
plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off')
panel_no +=1 
#if (i!=0):
#	plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off') 

mean = np.mean(stokes[:,:,0,:],axis=(0,1))
mean = flt.gaussian_filter(mean,20)
wls = argrelextrema(mean,np.less)
wls = np.asarray(wls)
wls = wls[0]
wls = np.append(0,wls)
print wls
ll = wls[-1]


#plt.subplot(panelsy,panelsx,2)
#plt.imshow(stokes[:,:,0,ll],origin='lower',cmap=Imap,vmin=0.2,vmax=0.6)
#if (i==2):
#plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
#plt.title('Stokes $I$ (line)')
#plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off') 
#plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off') 

plt.subplot(panelsy,panelsx,panel_no)
plt.imshow(np.mean(stokes[:,:,3,ll+4:ll+10],axis=2),origin='lower',cmap=Bmap,vmin=-0.05,vmax=0.05)
#if (i==2):
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Stokes $V$ (Iron line)')
plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off') 
plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off') 
panel_no +=1


# We will keep this simple:
for i in T_nodes:
	m = np.mean(parameters[i])
	s = np.std(parameters[i])
	k=k+1
	plt.subplot(panelsy,panelsx,panel_no)
	plt.imshow(parameters[i],origin='lower',cmap=Tmap,vmin=m-3*s,vmax=m+3*s)
	#if (i==2):
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.title('Temperature [K] at $\log\,\\tau$ = '+str(T_nodes_tau[i]))
	plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off') 
	if (i!=T_nodes[0]):
		plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off') 
	panel_no += 1


for i in vs_nodes:

	parameters[i] /= -1E5
	m = np.mean(parameters[i])
	s = np.std(parameters[i])

	plt.subplot(panelsy,panelsx,panel_no)
	plt.imshow(parameters[i],origin='lower',cmap=Vmap,vmin=-4*s,vmax=4*s)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.title('LOS velocity [km/s] at $\log\,\\tau$ = '+str(vs_nodes_tau[i-vs_nodes[0]]))
	plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off') 
	if (i!=vs_nodes[0]):
		plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off') 
	panel_no += 1


#s = np.copy(B_nodes)
for i in B_nodes:
	parameters[i] *= np.cos(parameters[theta_nodes[0]])
#	s[i-B_nodes[0]] = 3.0*np.std(parameters[i])

#for i in range(0,-1):
for i in B_nodes:
	print i
	plt.subplot(panelsy,panelsx,panel_no)
	plt.imshow(-parameters[i],origin='lower',cmap=Bmap,vmin=-1000,vmax=1000)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.title('$\mathrm{B\,[Gauss]}$ at $\log\,\\tau =$'+str(B_nodes_tau[i-B_nodes[0]]))
	if (i!=B_nodes[0]):
		plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off') 
	panel_no += 1


plt.tight_layout()
plt.savefig(print_maps_here+'.eps',fmt='eps',bbox_inches='tight')
plt.savefig(print_maps_here+'.png',fmt='png',bbox_inches='tight')
plt.clf()
plt.cla()

	






