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
cube1[-(NB+NT):-NT] *= cube1[-1]
cube2[-1] = np.cos(cube2[-1])
cube2[-(NB+NT):-NT] *= cube2[-1]

plt.clf()
plt.cla()

params = [0,1,2,3,4,5,6,7,8,9,10]
#cmaps = ['inferno','inferno','inferno','inferno','inferno','inferno','bwr','bwr','bwr','bwr','PRGn','PRGn']
#unit = ['K','K','km/s','km/s','Gauss','Gauss']
cmaps = [None] * 12
unit = [None] * 12
for i in range(0,4):
	cmaps[i] = 'inferno'
	unit[i] = 'K'
cmaps[4]  = 'inferno'
unit[4]  = 'km/s'
for i in range(5,8):
	cmaps[i] = 'bwr'
	unit[i] = 'km/s'
for i in range(8,11):
	cmaps[i] = 'PRGn'
	unit[i] = 'Gauss'
print cmaps

cubename1 = 'Low starting value'
cubename2 = 'High starting value'

index = 0
for p in params:
	if (cmaps[index] == 'bwr'):
		cube1[p] /= -1E5
		cube2[p] /= -1E5
	plt.figure(figsize=[12.5,3.5])
	plt.clf()
	plt.cla()
	plt.subplot(1,3,1)
	m = np.mean(cube1[p])
	s = np.std(cube1[p])
	plt.imshow(cube1[p],vmin=m-3*s,vmax=m+3*s,cmap=cmaps[index])
	plt.title(cubename1)
	plt.colorbar()
	plt.subplot(1,3,2)
	plt.imshow(cube2[p],vmin=m-3*s,vmax=m+3*s,cmap=cmaps[index])
	plt.title(cubename2)
	plt.colorbar()
	plt.subplot(1,3,3)
	#plt.hist(cube1[p].flatten(),bins=100,alpha=0.5,color='red')
	#plt.hist(cube2[p].flatten(),bins=100,alpha=0.5,color='blue')
	#plt.subplot(224)
	plt.hist2d(cube1[p].flatten(),cube2[p].flatten(),bins=(100,100),cmap='Purples',vmax=30,range=[[m-3*s,m+3*s],[m-3*s,m+3*s]])
	plt.plot(cube1[p].flatten(),cube1[p].flatten(),color='red')
	plt.xlabel(cubename1)
	plt.ylabel(cubename2)
	r = stats.pearsonr(cube1[p].flatten(),cube2[p].flatten())
	r = np.asarray(r)
	r[1] = np.std(cube1[p].flatten()-cube2[p].flatten())
	r[0] = round(r[0],3)
	r[1] = round(r[1],3)
	plt.text(m-3*s,m-2.8*s,'r='+str(r[0])+' $\\sigma=$'+str(r[1])+unit[index],fontsize=14)
	print r
	plt.tight_layout()
	plt.savefig(sys.argv[3]+'_'+str(p)+'_.png',fmt='png',bbox_inches='tight')
	plt.savefig(sys.argv[3]+'_'+str(p)+'_.eps',fmt='eps',bbox_inches='tight')
	index +=1

#Here we plot correlations of correlations

combo1 = [3,8]
combo2 = [8,11]
combo3 = [7,10]

plt.clf()
plt.cla()
plt.subplot(131)
v1 = (cube1[combo1[0]] * cube1[combo1[1]]).flatten()
v2 = (cube2[combo1[0]] * cube2[combo1[1]]).flatten()
m = np.mean(v1)
s = np.std(v1)
plt.hist2d(v1,v2,bins=(100,100),cmap='Purples',vmax=30,range=[[m-3*s,m+3*s],[m-3*s,m+3*s]])
plt.plot(v1,v1,color='red')
plt.xlabel('Standard inversion')
plt.ylabel('CNN')
r = stats.pearsonr(v1,v2)
plt.text(m-3*s,m-2.8*s,'r='+str(r[0]),fontsize=14)
plt.subplot(132)
v1 = (cube1[combo2[0]] * cube1[combo2[1]]).flatten()
v2 = (cube2[combo2[0]] * cube2[combo2[1]]).flatten()
m = np.mean(v1)
s = np.std(v1)
plt.hist2d(v1,v2,bins=(100,100),cmap='Purples',vmax=30,range=[[m-3*s,m+3*s],[m-3*s,m+3*s]])
plt.plot(v1,v1,color='red')
plt.xlabel('Standard inversion')
plt.ylabel('CNN')
r = stats.pearsonr(v1,v2)
plt.text(m-3*s,m-2.8*s,'r='+str(r[0]),fontsize=14)
plt.subplot(133)
v1 = (cube1[combo3[0]] * cube1[combo3[1]]).flatten()
v2 = (cube2[combo3[0]] * cube2[combo3[1]]).flatten()
m = np.mean(v1)
s = np.std(v1)
plt.hist2d(v1,v2,bins=(100,100),cmap='Purples',vmax=30,range=[[m-3*s,m+3*s],[m-3*s,m+3*s]])
plt.plot(v1,v1,color='red')
plt.xlabel('Standard inversion')
plt.ylabel('CNN')
r = stats.pearsonr(v1,v2)
plt.text(m-3*s,m-2.8*s,'r='+str(r[0]),fontsize=14)
plt.tight_layout()
plt.savefig('corr_of_corr.png',fmt='png',bbox_inches='tight')
plt.savefig('corr_of_corr.eps',fmt='eps',bbox_inches='tight')


	
	