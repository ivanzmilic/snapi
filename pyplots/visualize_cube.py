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

print_maps_here = sys.argv[3]

#Here we read the parameter/node map:
temp = pyana.fzread(input_nodes)
parameters = temp["data"]
temp = parameters.shape
NN = temp[0] #total number of nodes

a_read = pyana.fzread(input_atmos)
atmospheres = a_read["data"]


#-----------------------------------------------------------------------------------------------------
# NOW PLOT THE NODES ---------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------

T_nodes = [0,2,3]
T_nodes_tau = [-2.7,-1.8,-0.9,0.0]
#vt_nodes = [5]
vs_nodes = [4,5,6]
vs_nodes_tau = [-3.5,-1.7,-0.5]
B_nodes = [7,8,9]
B_nodes_tau = [-3.0,-1.7,-0.5]
theta_nodes = [10]

panelsx=3
panelsy=3

Tmap = 'hot'
Vmap = 'coolwarm'
Bmap = 'summer'
Imap = 'gray'
Pmap = 'Spectral'
Dmap = 'coolwarm'

plt.clf()
plt.cla()

defsize = 2.5
barshrink=0.8
plt.figure(figsize=[defsize*3*panelsx, defsize*panelsy])
k =0

# We will keep this simple:
for i in T_nodes:
	m = np.mean(parameters[i])
	s = np.std(parameters[i])
	k=k+1
	plt.subplot(panelsy,panelsx,k)
	plt.imshow(parameters[i],origin='lower',cmap=Tmap,vmin=m-3*s,vmax=m+3*s)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.title('Temperature at $\log\,\\tau$ = '+str(T_nodes_tau[i-T_nodes[0]]))


for i in range(vs_nodes[0],vs_nodes[-1]+1):

	parameters[i] /= -1E5
	m = np.mean(parameters[i])
	s = np.std(parameters[i])

	plt.subplot(panelsy,panelsx,panelsx+(i-vs_nodes[0])+1)
	plt.imshow(parameters[i],origin='lower',cmap=Vmap,vmin=-3*s,vmax=3*s)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.title('LOS velocity at $\log\,\\tau$ = '+str(vs_nodes_tau[i-vs_nodes[0]]))

s = np.copy(B_nodes)
for i in range(B_nodes[0],B_nodes[-1]+1):
	parameters[i] *= np.cos(parameters[theta_nodes[0]])
	s[i-B_nodes[0]] = 3.0*np.std(parameters[i])


for i in range(B_nodes[0],B_nodes[-1]+1):
	plt.subplot(panelsy,panelsx,2*panelsx+(i-B_nodes[0])+1)
	plt.imshow(-parameters[i],origin='lower',cmap=Bmap,vmin=0,vmax=1500)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.title('$\mathrm{B\,[Gauss]}$ at $\log\,\\tau =$'+str(B_nodes_tau[i-B_nodes[0]]))

plt.tight_layout()
plt.savefig(print_maps_here+'.eps',fmt='eps',bbox_inches='tight')
plt.savefig(print_maps_here+'.png',fmt='png',bbox_inches='tight')
plt.clf()
plt.cla()

	






