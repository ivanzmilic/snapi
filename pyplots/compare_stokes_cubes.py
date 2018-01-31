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
out_name = sys.argv[3]
maskfile = sys.argv[4]

temp = pyana.fzread(cube1_in)
cube1 = temp["data"]
#cube1 = np.transpose(cube1,(1,0,2,3))
temp = pyana.fzread(cube2_in)
cube2 = temp["data"]

l_l = 1
l_r = 501


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
plt.plot(cube_1_mean,color='red')
plt.plot(cube_2_mean,color='blue')
plt.savefig('mean_profiles',fmt='png')
cube_1_mean = flt.gaussian_filter(cube_1_mean,2)
wls = argrelextrema(cube_1_mean,np.less)
wls = np.asarray(wls)
wls = wls[0]
wls = np.append(0,wls)
print wls
#wls[0] += 15
#wls[3] += 15
N_x_panels = 6
N_y_panels = len(wls)

#make the size of the figure:
x_size = 4.0
y_size = x_size * float(NY)/float(NX) * 1.0

shrinkage = 0.6

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

	to_plot1 = np.mean(cube1[:,:,3,wls[j-1]:wls[j-1]+6],axis=2)/cube1[:,:,0,0]
	s = np.std(to_plot1)

	plt.subplot(N_y_panels,N_x_panels,(j-1)*N_x_panels+4)
	plt.imshow(to_plot1,origin='lower',vmin = -3*s,vmax=3*s,cmap='coolwarm')
	plt.colorbar(shrink=shrinkage)

	to_plot2 = np.mean(cube2[:,:,3,wls[j-1]:wls[j-1]+6],axis=2)/cube1[:,:,0,0]
	
	plt.subplot(N_y_panels,N_x_panels,(j-1)*N_x_panels+5)
	plt.imshow(to_plot2,origin='lower',vmin = -3*s,vmax=3*s,cmap='coolwarm')
	plt.colorbar(shrink=shrinkage)

	to_plot = (to_plot2-to_plot1)/(to_plot1+1E-3*np.amax(to_plot1))*100.0
	plt.subplot(N_y_panels,N_x_panels,(j-1)*N_x_panels+6)
	plt.imshow(to_plot,origin='lower',vmin = -20.0,vmax=20.0,cmap='coolwarm')
	plt.colorbar(shrink=shrinkage)


plt.tight_layout()
plt.savefig(out_name,fmt='png',bbxox_inches='tight')
plt.savefig(out_name+'.eps',fmt='eps',bbxox_inches='tight')
