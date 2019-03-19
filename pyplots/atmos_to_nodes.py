# a small piece of code, which interpolates input atmosphere
# to a very sparse space of nodes in order to get a node-based model from
# given model atmosphere. Used when we want to change node configuration 
# during the fitting.

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import pyana
import numpy as np 
import sys
import scipy.interpolate as interpol

#read the atmosphere, assumes f0 format:
atmos_in = sys.argv[1]
temp = pyana.fzread(atmos_in)
atmos = temp["data"]

# now we need the nodes:
T_nodes = [-3.0,-2.1,-1.3,-6.0,0.0]
vt_nodes = [0]
v_nodes = [-3.3,-0.5]
B_nodes = [-2.8,0.3]
theta_nodes=[0]
phi_nodes=[]

#now time to interpolate
NX = atmos.shape[0]
NY = atmos.shape[1]

N_T = len(T_nodes)
N_vt = len(vt_nodes)
N_v = len(v_nodes)
N_B = len(B_nodes)
N_theta = len(theta_nodes)
N_phi = len(phi_nodes)
N_nodes = N_T + N_vt + N_v + N_B + N_theta + N_phi
nodes = np.zeros([N_nodes,NX,NY])


for i in range(0,NX):
	for j in range(0,NY):

		index = 0
		f = interpol.interp1d(atmos[i,j,0],atmos[i,j,2])
		T = f(T_nodes)
		nodes[index:index+N_T,i,j] = T
		index+=N_T
		f = interpol.interp1d(atmos[i,j,0],atmos[i,j,8])
		vt = f(vt_nodes)
		nodes[index:index+N_vt,i,j] = vt
		index+=N_vt
		f = interpol.interp1d(atmos[i,j,0],atmos[i,j,9])
		v = f(v_nodes)
		nodes[index:index+N_v,i,j] = v
		index+=N_v
		f = interpol.interp1d(atmos[i,j,0],atmos[i,j,7])
		B = f(B_nodes)
		nodes[index:index+N_B,i,j] = B
		index+=N_B
		f = interpol.interp1d(atmos[i,j,0],atmos[i,j,10])
		theta = f(theta_nodes)
		nodes[index:index+N_theta,i,j] = theta
		index+=N_theta
		f = interpol.interp1d(atmos[i,j,0],atmos[i,j,11])
		phi = f(phi_nodes)
		nodes[index:index+N_phi,i,j] = phi

#additional polishing - hardcode
#nodes[4,:,:] = 2E4

nodes = np.transpose(nodes,(0,2,1))
#write down nodes in the file:
output_file = sys.argv[2]
pyana.fzwrite(output_file,nodes,0,'placeholder')
		


