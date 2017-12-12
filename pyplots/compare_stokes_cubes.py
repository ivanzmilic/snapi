import matplotlib
matplotlib.use('Agg')
import pyana
import numpy as np 
import matplotlib.pyplot as plt 
import sys
from scipy.signal import argrelextrema

cube1_in = sys.argv[1]
cube2_in = sys.argv[2]
out_name = sys.argv[3]
maskfile = sys.argv[4]

temp = pyana.fzread(cube1_in)
cube1 = temp["data"]
cube1 = np.transpose(cube1,(1,0,2,3))
temp = pyana.fzread(cube2_in)
cube2 = temp["data"]

l_l = 1
l_r = 601


dims = cube1.shape
NX = dims[1]
NY = dims[0]
NL = dims[3]
print NX, NY

#before plotting the spatial distributions, plot averaged spectra:

cube_1_mean = np.mean(cube1[:,:,0,:],axis=(0,1))
cube_2_mean = np.mean(cube2[:,:,0,:],axis=(0,1))
plt.clf()
plt.cla()
plt.plot(cube_1_mean)
plt.plot(cube_2_mean)
plt.savefig('mean_profiles',fmt='png')
wls = argrelextrema(cube_2_mean,np.less)
wls = np.asarray(wls)
wls = wls[0]
print wls
#wls[0] += 15
#wls[3] += 15
N_x_panels = 3
N_y_panels = len(wls)

#make the size of the figure:
x_size = 6.0
y_size = x_size * float(NY)/float(NX) * 1.0

shrinkage = 1.0

plt.figure(figsize=[N_x_panels*x_size,N_y_panels*y_size])

for j in range (1,N_y_panels+1):
	
	to_plot = np.copy(cube1[:,:,0,wls[j-1]])
	m = np.mean(to_plot)
	to_plot/=m
	s = np.std(to_plot)
	
	plt.subplot(N_y_panels,N_x_panels,(j-1)*N_x_panels+1)
	plt.imshow(to_plot,origin='lower',vmin = 1-3*s,vmax=1+3*s,cmap='hot')
	plt.colorbar(shrink=shrinkage)

	to_plot = np.copy(cube2[:,:,0,wls[j-1]])
	to_plot/=m
	plt.subplot(N_y_panels,N_x_panels,(j-1)*N_x_panels+2)
	plt.imshow(to_plot,origin='lower',vmin = 1-3*s,vmax=1+3*s,cmap='hot')
	plt.colorbar(shrink=shrinkage)


	to_plot = (cube2[:,:,0,wls[j-1]]-cube1[:,:,0,wls[j-1]])/cube1[:,:,0,wls[j-1]]*100.0
	plt.subplot(N_y_panels,N_x_panels,(j-1)*N_x_panels+3)
	plt.imshow(to_plot,origin='lower',vmin = -20.0,vmax=20.0,cmap='coolwarm')
	plt.colorbar(shrink=shrinkage)

plt.tight_layout()
plt.savefig(out_name,fmt='png',bbxox_inches='tight')
plt.savefig(out_name+'.eps',fmt='eps',bbxox_inches='tight')
