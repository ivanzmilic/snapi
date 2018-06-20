import matplotlib
matplotlib.use('Agg')
import pyana
import numpy as np 
import matplotlib.pyplot as plt 
import sys
from scipy.signal import argrelextrema
import scipy.ndimage.filters as flt

cube1_in = sys.argv[1]
cube2_in = sys.argv[2]
nodes_in = sys.argv[3]

temp = pyana.fzread(cube1_in)
cube1 = temp["data"]
temp = pyana.fzread(cube2_in)
cube2 = temp["data"]

temp = pyana.fzread(nodes_in)
nodes = temp["data"]

l_l = 1
l_r = 501


dims = cube1.shape
NX = dims[1]
NY = dims[0]
NL = dims[3]
print NX, NY

chisq = (cube1-cube2)
chisq = chisq*chisq
chisq = np.sum(chisq[:,:,0,:],axis=2) + np.sum(chisq[:,:,3,:],axis=2)
chisq /= 1001 * 1E25

chi_hist = chisq.reshape(NX*NY)

plt.cla()
plt.clf()
plt.hist(chi_hist,bins=50)
plt.savefig('chisq_hist',fmt='png')

plt.cla()
plt.clf()
plt.imshow(chisq,origin='lower',vmin=0,vmax=10)
plt.colorbar()
plt.savefig('chisq_map',fmt='png')

for i in range(0,NX):
	for j in range(0,NY):
		chisq_min = chisq[i,j]
		if chisq[i,j] > 1.0:
			for ii in range(0,NX):
				for jj in range(0,NY):
					res = cube1[i,j]-cube2[ii,jj]
					res = res*res
					chisq_loc = np.sum(res[0]) + np.sum(res[3])
					chisq_loc /= 1001 * 1E25
					if chisq_loc < chisq_min:
						chisq_min = chisq_loc
						nodes[:,i,j] = np.copy(nodes[:,ii,jj])
						cube2[i,j] = np.copy(cube2[ii,jj])
						
		chisq[i,j] = chisq_min
	print 'i = ',i,'done'

chi_hist = chisq.reshape(NX*NY)
plt.cla()
plt.clf()
plt.hist(chi_hist,bins=50)
plt.savefig('chisq_hist_new',fmt='png')
pyana.fzwrite('inverted_spectra_purified.f0',cube2,0,'placeholder')
pyana.fzwrite('nodes_purified.f0',nodes,0,'placeholder')


