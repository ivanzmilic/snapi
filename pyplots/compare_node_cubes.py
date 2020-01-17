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

print ('Number of parameters : ', NP)

NB = 2
NT = 1
cube1[-1] = np.cos(cube1[-1])
cube1[-(NB+NT):-NT] *= cube1[-1]
cube2[-1] = np.cos(cube2[-1])
cube2[-(NB+NT):-NT] *= cube2[-1]

plt.clf()
plt.cla()

params = [0,1,2,3,4,5,6,7,8,9,10,11]
#cmaps = ['inferno','inferno','inferno','inferno','inferno','inferno','bwr','bwr','bwr','bwr','PRGn','PRGn']
#unit = ['K','K','km/s','km/s','Gauss','Gauss']
cmaps = [None] * 13
unit = [None] * 13
for i in range(0,5):
	cube1[i] /= 1E3
	cube2[i] /= 1E3
	cmaps[i] = 'inferno'
	unit[i] = 'kK'
cmaps[5]  = 'inferno'
unit[5]  = 'km/s'
for i in range(6,10):
	cmaps[i] = 'bwr'
	unit[i] = 'km/s'
for i in range(10,12):
	cmaps[i] = 'PRGn'
	unit[i] = 'Gauss'
print (cmaps)

cubename1 = 'Standard inversion'
cubename2 = 'CNN'

x = np.linspace(0,287,288)*20.8/1E3
y = np.linspace(0,287,288)*20.8/1E3

index = 0
for p in params:
	if (cmaps[index] == 'bwr'):
		cube1[p] /= -1E5
		cube2[p] /= -1E5
	plt.figure(figsize=[12.0,4.0])
	plt.clf()
	plt.cla()
	plt.subplot(1,3,1)
	m = np.mean(cube1[p])
	s = np.std(cube1[p])
	plt.imshow(cube1[p],vmin=m-3*s,vmax=m+3*s,cmap=cmaps[index],origin='lower',extent=[x[0],x[-1],y[0],y[-1]])
	plt.title(cubename1+' ['+str(unit[p])+']')
	plt.xlabel('$x\,\mathrm{[Mm]}$')
	plt.ylabel('$y\,\mathrm{[Mm]}$')
	plt.colorbar(shrink=0.8)
	plt.subplot(1,3,2)
	plt.imshow(cube2[p],vmin=m-3*s,vmax=m+3*s,cmap=cmaps[index],origin='lower',extent=[x[0],x[-1],y[0],y[-1]])
	plt.xlabel('$x\,\mathrm{[Mm]}$')
	plt.title(cubename2+' ['+str(unit[p])+']')
	plt.colorbar(shrink=0.8)
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
	plt.text(m-3*s,m-2.8*s,'r='+str(r[0])+' $\\sigma=$'+str(r[1])+''+unit[index],fontsize=14)
	print (r)
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
plt.figure(figsize=[6.0,2.0])
plt.subplot(121)
v1 = (cube1[combo1[0]] * cube1[combo1[1]]).flatten()
v2 = (cube2[combo1[0]] * cube2[combo1[1]]).flatten()
m = np.mean(v1)
s = np.std(v1)
plt.hist2d(v1,v2,bins=(100,100),cmap='Purples',vmax=30,range=[[m-3*s,m+3*s],[m-3*s,m+3*s]])
plt.plot(v1,v1,color='red')
plt.xlabel('Standard inversion')
plt.ylabel('CNN')
r = stats.pearsonr(v1,v2)
r = np.asarray(r)
r[0] = round(r[0],3)
plt.text(m-3*s,m-2.8*s,'r='+str(r[0]),fontsize=14)
plt.subplot(122)
v1 = (cube1[combo2[0]] * cube1[combo2[1]]).flatten()
v2 = (cube2[combo2[0]] * cube2[combo2[1]]).flatten()
m = np.mean(v1)
s = np.std(v1)
plt.hist2d(v1,v2,bins=(100,100),cmap='Purples',vmax=30,range=[[m-3*s,m+3*s],[m-3*s,m+3*s]])
plt.plot(v1,v1,color='red')
plt.xlabel('Standard inversion')
r = stats.pearsonr(v1,v2)
r = np.asarray(r)
r[0] = round(r[0],3)
plt.text(m-3*s,m-2.8*s,'r='+str(r[0]),fontsize=14)
plt.tight_layout()
plt.savefig('corr_of_corr.png',fmt='png',bbox_inches='tight')
plt.savefig('corr_of_corr.eps',fmt='eps',bbox_inches='tight')


	
	