import numpy as np 
import scipy.signal as sgn
import matplotlib.pyplot as plt 
import sys 
import pyana 
import matplotlib

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
temp = parameters.shape
NN = temp[0] #total number of nodes

a_read = pyana.fzread(input_atmos)
atmospheres = a_read["data"]

temp = pyana.fzread(input_spectra)
stokes = temp["data"]

I_c = np.mean(stokes[:,:,0,0])
stokes[:,:,:,:] /= I_c

#-----------------------------------------------------------------------------------------------------
# NOW PLOT THE NODES ---------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------

T_nodes = [0,2,3]
T_nodes_tau = [-3.0,-1.7,-0.7,0.0]
#vt_nodes = [3]
vs_nodes = [4,5,6]
vs_nodes_tau = [-3.5,-1.7,-0.5]
B_nodes = [7,8,9]
B_nodes_tau = [-3.0,-1.6,-0.5]
theta_nodes = [10]

panelsx=3
panelsy=4

Tmap = 'hot'
Vmap = 'coolwarm'
Bmap = 'summer'
Imap = 'hot'
Pmap = 'Spectral'
Dmap = 'coolwarm'

plt.clf()
plt.cla()

defsize = 2.5
barshrink=0.8
plt.figure(figsize=[defsize*2.5*panelsx, 0.95*defsize*panelsy])
k =0

plt.subplot(panelsy,panelsx,1)
plt.imshow(stokes[:,:,0,0].transpose(),origin='lower',cmap=Imap,vmin=0.7,vmax=1.3)
#if (i==2):
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Stokes $I$ (continuum)')
plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off') 
#if (i!=0):
#	plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off') 

ll =240

plt.subplot(panelsy,panelsx,2)
plt.imshow(stokes[:,:,0,ll].transpose(),origin='lower',cmap=Imap,vmin=0.2,vmax=0.6)
#if (i==2):
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Stokes $I$ (Sodium D2)')
plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off') 
plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off') 

plt.subplot(panelsy,panelsx,3)
plt.imshow(stokes[:,:,3,ll].transpose(),origin='summer',cmap=Bmap,vmin=0,vmax=0.03)
#if (i==2):
plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
plt.title('Stokes $V$ (Sodium D2)')
plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off') 
plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off') 


# We will keep this simple:
for i in T_nodes:
	m = np.mean(parameters[i])
	s = np.std(parameters[i])
	k=k+1
	plt.subplot(panelsy,panelsx,panelsx+k)
	plt.imshow(parameters[i],origin='lower',cmap=Tmap,vmin=m-3*s,vmax=m+3*s)
	#if (i==2):
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.title('Temperature [K] at $\log\,\\tau$ = '+str(T_nodes_tau[i-T_nodes[0]]))
	plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off') 
	if (i!=0):
		plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off') 


for i in range(vs_nodes[0],vs_nodes[-1]+1):

	parameters[i] /= -1E5
	m = np.mean(parameters[i])
	s = np.std(parameters[i])

	plt.subplot(panelsy,panelsx,2*panelsx+(i-vs_nodes[0])+1)
	plt.imshow(parameters[i],origin='lower',cmap=Vmap,vmin=-4*s,vmax=4*s)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.title('LOS velocity [km/s] at $\log\,\\tau$ = '+str(vs_nodes_tau[i-vs_nodes[0]]))
	plt.tick_params(axis='x',which='both',bottom='off',top='off',labelbottom='off') 
	if (i!=4):
		plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off') 


s = np.copy(B_nodes)
for i in range(B_nodes[0],B_nodes[-1]+1):
	parameters[i] *= np.cos(parameters[theta_nodes[0]])
	s[i-B_nodes[0]] = 3.0*np.std(parameters[i])

#for i in range(0,-1):
for i in range(B_nodes[0],B_nodes[-1]+1):
	plt.subplot(panelsy,panelsx,3*panelsx+(i-B_nodes[0])+1)
	plt.imshow(-parameters[i],origin='lower',cmap=Bmap,vmin=0,vmax=1500)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.title('$\mathrm{B\,[Gauss]}$ at $\log\,\\tau =$'+str(B_nodes_tau[i-B_nodes[0]]))
	if (i!=7):
		plt.tick_params(axis='y',which='both',left='off',right='off',labelleft='off') 


plt.tight_layout()
plt.savefig(print_maps_here+'.eps',fmt='eps',bbox_inches='tight')
plt.savefig(print_maps_here+'.png',fmt='png',bbox_inches='tight')
plt.clf()
plt.cla()

	






