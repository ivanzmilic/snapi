import pyana
import matplotlib.pyplot as plt
import numpy as np 
import scipy.stats as stats

import sys

file1 = sys.argv[1]
file2 = sys.argv[2] 

cube1 = pyana.fzread(file1)["data"]
cube2 = pyana.fzread(file2)["data"]

NP = cube1.shape[0]

print 'Number of parameters : ', NP

NB = 3
NT = 1
cube1[-1] = np.cos(cube1[-1])
cube1[-NB-NT:-NT] *= cube1[-1]
cube2[-1] = np.cos(cube2[-1])
cube2[-NB-NT:-NT] *= cube2[-1]

plt.clf()
plt.cla()
plt.figure(figsize=[16.5,4*NP])

params = [3,4,8,9,10,11]
cmaps = ['magma','magma','bwr']

for p in range(0,NP):
	plt.subplot(NP,3,3*p+1)
	m = np.mean(cube1[p])
	s = np.std(cube1[p])
	plt.imshow(cube1[p],vmin=m-3*s,vmax=m+3*s)
	plt.colorbar()
	plt.subplot(NP,3,3*p+2)
	plt.imshow(cube2[p],vmin=m-3*s,vmax=m+3*s)
	plt.colorbar()
	plt.subplot(NP,3,3*p+3)
	#plt.hist(cube1[p].flatten(),bins=100,alpha=0.5,color='red')
	#plt.hist(cube2[p].flatten(),bins=100,alpha=0.5,color='blue')
	#plt.subplot(224)
	plt.hist2d(cube1[p].flatten(),cube2[p].flatten(),bins=(100,100),cmap='Purples',vmax=30,range=[[m-3*s,m+3*s],[m-3*s,m+3*s]])
	plt.plot(cube1[p].flatten(),cube1[p].flatten(),color='red')
	r = stats.pearsonr(cube1[p].flatten(),cube2[p].flatten())
	r = np.asarray(r)
	r[1] = np.std(cube1[p].flatten()-cube2[p].flatten())
	plt.text(m-3*s,m-3*s,str(r))
	print r

plt.tight_layout()
plt.savefig(sys.argv[3]+'.png',fmt='png',bbox_inches='tight')
	
	