import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyana
import numpy as np 
import sys
from scipy.interpolate import interp1d

node_in = sys.argv[1]
cube_in = sys.argv[2]
plot_here = sys.argv[3]

temp = pyana.fzread(node_in)
nodes = temp["data"]

dim = nodes.shape
NX = dim[1]
NY = dim[2]

temp = pyana.fzread(cube_in)
cube = temp["data"]

print nodes.shape
print cube.shape

plt.clf()
plt.cla()
plt.figure(figsize=[9,9])

c_map = 'coolwarm'
c_map = 'hot'


tau = [-1.5,-0.5,0.5]
tau = [-1.5,-0.4,0.3]
tau = np.asarray(tau)
n = [0,1,2]
n_nodes = len(tau)

cube_interpol = np.zeros([n_nodes,NX,NY])

p = 9
p = 2

for i in range(0,NX):
	for j in range(0,NY):
		f = interp1d(cube[i,j,0,:],cube[i,j,p,:])
		cube_interpol[:,i,j] = f(tau)

for i in range (0,n_nodes):

	s = np.std(cube_interpol[i,:,:])
	m = np.mean(cube_interpol[i,:,:])
	print s,m
	plt.subplot(n_nodes,3,i*3+1)
	plt.imshow(cube_interpol[i,:,:],origin='lower',cmap=c_map,vmin=m-3*s,vmax=m+3*s)
	plt.title('Simulation')
	
	plt.subplot(n_nodes,3,i*3+2)
	plt.imshow(nodes[n[i]],origin='lower',cmap=c_map,vmin=m-3*s,vmax=m+3*s)
	plt.title('Inversion')
	plt.ylabel('Node # '+str(n[i]))
	
	plt.subplot(n_nodes,3,i*3+3)
	#plt.imshow((nodes[n] - cube[:144,:144,9,d])/cube[:144,:144,9,d],origin='lower',vmin=-0.3,vmax=0.3,cmap=c_map)
	#diff = (nodes[n[i]] - cube[x_low:x_high,y_low:y_high,7,d[i]])
	#diff = diff.reshape((x_high-x_low)*(y_high-y_low))
	#plt.hist(diff,50,range=[-100,100], normed=1, facecolor='green', alpha=0.75)
	#plt.title('Histogram of differences')

plt.tight_layout()

plt.savefig(plot_here,fmt='eps',bbox_inches='tight')
plt.savefig(plot_here,fmt='png',bbox_inches='tight')








