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


to_plot = [2,3,8,10,11]
NP = len(to_plot)
plt.clf()
plt.cla()
scale = 3.4
ratio = 0.93
fig, axes = plt.subplots(nrows=NP,ncols=3,figsize=[scale*3,scale*NP*ratio])
image_no = 0

for i in range (0,len(to_plot)):
	
	p = to_plot[i]
	# if it is the los velocity, scale and orient it properly
	if (cmaps[p] == 'bwr'):
		cube1[p] /= -1E5
		cube2[p] /= -1E5
	
	m = np.mean(cube1[p])
	s = np.std(cube1[p])
	
	ax = axes.flat[image_no]
	im = ax.imshow(cube1[p],vmin=m-3*s,vmax=m+3*s,cmap=cmaps[p],origin='lower',extent=[x[0],x[-1],y[0],y[-1]])
	if (i==0):
		ax.set_title('Standard Inversion')
	if (i==NP-1):
		ax.set_xlabel('$x\,\mathrm{[Mm]}$')
	ax.set_ylabel('$y\,\mathrm{[Mm]}$')
	if (i!=NP-1):
		ax.set_xticklabels([])
	
	image_no +=1
	ax = axes.flat[image_no]
	im = ax.imshow(cube2[p],vmin=m-3*s,vmax=m+3*s,cmap=cmaps[p],origin='lower',extent=[x[0],x[-1],y[0],y[-1]])
	if (i==0):
		ax.set_title('CNN')
	if (i==NP-1):
		ax.set_xlabel('$x\,\mathrm{[Mm]}$')
	if (i!=NP-1):
		ax.set_xticklabels([])
	ax.set_yticklabels([])
	ax.colorbar(shrink=0.8)

	image_no +=1
	ax = axes.flat[image_no]
	im = ax.hist2d(cube1[p].flatten(),cube2[p].flatten(),bins=(100,100),cmap='Purples',vmax=30,range=[[m-3*s,m+3*s],[m-3*s,m+3*s]])
	ax.plot(cube1[p].flatten(),cube1[p].flatten(),color='red')
	if (i==NP-1):
		ax.set_xlabel(cubename1)
	
	ax.set_ylabel(cubename2)
	image_no +=1
	#r = stats.pearsonr(cube1[p].flatten(),cube2[p].flatten())
	#r = np.asarray(r)
	difference = (cube1[p].flatten()-cube2[p].flatten())
	median = np.percentile(difference,50)
	up = np.percentile(difference,95) - median
	down = median - np.percentile(difference,5)
	print (p,median,down,up)
	#r[0] = round(r[0],3)
	#r[1] = round(r[1],3)
	#plt.text(m-3*s,m-2.8*s,'r='+str(r[0])+' $\\sigma=$'+str(r[1])+''+unit[index],fontsize=14)
	#print (r)
fig.tight_layout()
fig.savefig(sys.argv[3]+'.png',fmt='png',bbox_inches='tight')
fig.savefig(sys.argv[3]+'.eps',fmt='eps',bbox_inches='tight')

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


	
	