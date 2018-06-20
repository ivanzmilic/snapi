import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyana
import numpy as np 
import sys

node_in = sys.argv[1]
cube_in = sys.argv[2]
plot_here = sys.argv[3]

temp = pyana.fzread(node_in)
nodes = temp["data"]

temp = pyana.fzread(cube_in)
cube = temp["data"]

print nodes.shape
print cube.shape

plt.clf()
plt.cla()
plt.figure(figsize=[9,9])

c_map = 'coolwarm'

y_low = 125
y_high = 200
x_low = 0
x_high = 75

tau = [-3.5,-1.6,-0.5]
tau = np.asarray(tau)
n = [8,9,10]
d = 70 + tau/0.1
d = np.asarray(d)
d.astype(int)
print d
#d = [41,56,65]

n_nodes = 3

nodes[4:8] /= 1E5
nodes[8:11] *= np.cos(nodes[11])
cube[:,:,7,:] *= np.cos(cube[:,:,10,:])

for i in range (0,n_nodes):

	plt.subplot(n_nodes,3,i*3+1)
	plt.imshow(nodes[n[i]],origin='lower',cmap=c_map,vmin=-1000,vmax=1000)
	plt.title('Inversion')
	plt.ylabel('Node # '+str(n[i]))
	print np.mean(nodes[n[i]])
	
	plt.subplot(n_nodes,3,i*3+2)
	plt.imshow(cube[x_low:x_high,y_low:y_high,7,d[i]],origin='lower',cmap=c_map,vmin=-1000,vmax=1000)
	plt.title('Simulation')
	print np.mean(cube[x_low:x_high,y_low:y_high,7,d[i]])
	
	plt.subplot(n_nodes,3,i*3+3)
	#plt.imshow((nodes[n] - cube[:144,:144,9,d])/cube[:144,:144,9,d],origin='lower',vmin=-0.3,vmax=0.3,cmap=c_map)
	diff = (nodes[n[i]] - cube[x_low:x_high,y_low:y_high,7,d[i]])
	diff = diff.reshape((x_high-x_low)*(y_high-y_low))
	plt.hist(diff,50,range=[-100,100], normed=1, facecolor='green', alpha=0.75)
	plt.title('Histogram of differences')

plt.tight_layout()

plt.savefig(plot_here,fmt='eps',bbox_inches='tight')
plt.savefig(plot_here,fmt='png',bbox_inches='tight')








