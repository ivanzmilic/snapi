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
#c_map = 'hot'


tau = [-1.5,-0.5,0.5]
#tau = [-1.5,-0.4,0.3]
tau = np.asarray(tau)
n = [5,6,7]
n_nodes = len(tau)

cube_interpol = np.zeros([n_nodes,NX,NY])
sign = -1

p = 9
#p = 2

for i in range(0,NX):
	for j in range(0,NY):
		f = interp1d(cube[i,j,0,:],cube[i,j,p,:])
		cube_interpol[:,i,j] = f(tau)

for i in range (0,n_nodes):

	s = np.std(cube_interpol[i,:,:])
	m = np.mean(cube_interpol[i,:,:])
	print s,m
	plt.subplot(n_nodes,2,i*2+1)
	plt.imshow(cube_interpol[i,:,:]*sign,origin='lower',cmap=c_map,vmin=m-3*s,vmax=m+3*s)
	plt.title('Simulation')
	
	plt.subplot(n_nodes,2,i*2+2)
	plt.imshow(nodes[n[i]]*sign,origin='lower',cmap=c_map,vmin=m-3*s,vmax=m+3*s)
	plt.title('Inversion')
	plt.ylabel('Node # '+str(n[i]))

	print np.std(nodes[n[i]]*sign - cube_interpol[i,:,:]*sign) / np.abs(np.amax(cube_interpol[i,:,:]*sign))

plt.tight_layout()

plt.savefig(plot_here+'_V',fmt='eps',bbox_inches='tight')
plt.savefig(plot_here+'_V',fmt='png',bbox_inches='tight')


tau = [-2.0,-0.9,-0.0,0.5]
n = [0,1,2,3]
n_nodes = len(tau)
cube_interpol = np.zeros([n_nodes,NX,NY])
sign = 1
p = 2
for i in range(0,NX):
	for j in range(0,NY):
		f = interp1d(cube[i,j,0,:],cube[i,j,p,:])
		cube_interpol[:,i,j] = f(tau)

plt.clf()
plt.cla()
plt.figure(figsize=[9,9])

c_map = 'hot'
for i in range (0,n_nodes):

	s = np.std(cube_interpol[i,:,:])
	m = np.mean(cube_interpol[i,:,:])
	print s,m
	plt.subplot(n_nodes,2,i*2+1)
	plt.imshow(cube_interpol[i,:,:]*sign,origin='lower',cmap=c_map,vmin=m-3*s,vmax=m+3*s)
	plt.title('Simulation')
	
	plt.subplot(n_nodes,2,i*2+2)
	plt.imshow(nodes[n[i]]*sign,origin='lower',cmap=c_map,vmin=m-3*s,vmax=m+3*s)
	plt.title('Inversion')
	plt.ylabel('Node # '+str(n[i]))
	
	#plt.subplot(n_nodes,3,i*3+3)
	#plt.imshow((nodes[n] - cube[:144,:144,9,d])/cube[:144,:144,9,d],origin='lower',vmin=-0.3,vmax=0.3,cmap=c_map)
	#diff = (nodes[n[i]] - cube[x_low:x_high,y_low:y_high,7,d[i]])
	#diff = diff.reshape((x_high-x_low)*(y_high-y_low))
	#plt.hist(diff,50,range=[-100,100], normed=1, facecolor='green', alpha=0.75)
	#plt.title('Histogram of differences')
	print np.std(nodes[n[i]]*sign - cube_interpol[i,:,:]*sign) / np.abs(np.amax(cube_interpol[i,:,:]*sign))

plt.tight_layout()

plt.savefig(plot_here+'_T',fmt='eps',bbox_inches='tight')
plt.savefig(plot_here+'_T',fmt='png',bbox_inches='tight')

tau = [-1.5,0.3]
n = [8,9]
n_nodes = len(tau)
cube_interpol = np.zeros([n_nodes,NX,NY])
sign = 1
p = 7
cube[:,:,7,:] *= np.cos(cube[:,:,10,:])
for i in range(0,NX):
	for j in range(0,NY):
		f = interp1d(cube[i,j,0,:],cube[i,j,p,:])
		cube_interpol[:,i,j] = f(tau)

plt.clf()
plt.cla()
plt.figure(figsize=[9,9])

nodes[n] *= np.cos(nodes[-1])

c_map = 'coolwarm'
for i in range (0,n_nodes):

	s = np.std(cube_interpol[i,:,:])
	m = np.mean(cube_interpol[i,:,:])
	print s,m
	plt.subplot(n_nodes,2,i*2+1)
	plt.imshow(cube_interpol[i,:,:]*sign,origin='lower',cmap=c_map,vmin=m-3*s,vmax=m+3*s)
	plt.title('Simulation')
	
	plt.subplot(n_nodes,2,i*2+2)
	plt.imshow(nodes[n[i]]*sign,origin='lower',cmap=c_map,vmin=m-3*s,vmax=m+3*s)
	plt.title('Inversion')
	plt.ylabel('Node # '+str(n[i]))
	
	#plt.subplot(n_nodes,3,i*3+3)
	#plt.imshow((nodes[n] - cube[:144,:144,9,d])/cube[:144,:144,9,d],origin='lower',vmin=-0.3,vmax=0.3,cmap=c_map)
	#diff = (nodes[n[i]] - cube[x_low:x_high,y_low:y_high,7,d[i]])
	#diff = diff.reshape((x_high-x_low)*(y_high-y_low))
	#plt.hist(diff,50,range=[-100,100], normed=1, facecolor='green', alpha=0.75)
	#plt.title('Histogram of differences')
	print np.std(nodes[n[i]]*sign - cube_interpol[i,:,:]*sign) / np.abs(np.amax(cube_interpol[i,:,:]*sign))

plt.tight_layout()

plt.savefig(plot_here+'_B',fmt='eps',bbox_inches='tight')
plt.savefig(plot_here+'_B',fmt='png',bbox_inches='tight')








