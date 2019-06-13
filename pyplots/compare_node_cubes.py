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

params = [2,3,7,8,10,11]
cmaps = ['inferno','inferno','bwr','bwr','PRGn','PRGn']

index = 0
for p in params:
	if (cmaps[index] == 'bwr'):
		cube1[p] /= -1E5
		cube2[p] /= -1E5
	plt.figure(figsize=[13.0,3.5])
	plt.clf()
	plt.cla()
	plt.subplot(1,3,1)
	m = np.mean(cube1[p])
	s = np.std(cube1[p])
	plt.imshow(cube1[p],vmin=m-3*s,vmax=m+3*s,cmap=cmaps[index])
	plt.title('Standard Inversion')
	plt.colorbar()
	plt.subplot(1,3,2)
	plt.imshow(cube2[p],vmin=m-3*s,vmax=m+3*s,cmap=cmaps[index])
	plt.title('CNN')
	plt.colorbar()
	plt.subplot(1,3,3)
	#plt.hist(cube1[p].flatten(),bins=100,alpha=0.5,color='red')
	#plt.hist(cube2[p].flatten(),bins=100,alpha=0.5,color='blue')
	#plt.subplot(224)
	plt.hist2d(cube1[p].flatten(),cube2[p].flatten(),bins=(100,100),cmap='Purples',vmax=30,range=[[m-3*s,m+3*s],[m-3*s,m+3*s]])
	plt.plot(cube1[p].flatten(),cube1[p].flatten(),color='red')
	r = stats.pearsonr(cube1[p].flatten(),cube2[p].flatten())
	r = np.asarray(r)
	r[1] = np.std(cube1[p].flatten()-cube2[p].flatten())
	plt.text(m-3*s,m-3*s,str(r),fontsize=18)
	print r
	plt.tight_layout()
	plt.savefig(sys.argv[3]+'_'+str(p)+'_.png',fmt='png',bbox_inches='tight')
	index +=1


	
	