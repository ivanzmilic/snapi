from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

import matplotlib
matplotlib.use('Agg')
import pyana
import numpy as np 
import matplotlib.pyplot as plt 
import sys
from scipy.signal import argrelextrema
import scipy.ndimage.filters as flt
from matplotlib import ticker
import colorcet as cc

cube1_in = sys.argv[1]
cube2_in = sys.argv[2]
out_name = sys.argv[3]
ifmask = sys.argv[4]
maskfile = sys.argv[5]

temp = pyana.fzread(cube1_in)
cube1 = temp["data"]
#cube1 = np.transpose(cube1,(1,0,2,3))
temp = pyana.fzread(cube2_in)
cube2 = temp["data"]

if (int(ifmask)):
	mask = np.loadtxt(maskfile,skiprows=1)
	print mask.shape


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
if (int(ifmask)):
	plt.plot(mask*cube_1_mean[0],'*')
plt.savefig('mean_profiles',fmt='png')
cube_1_mean = flt.gaussian_filter(cube_1_mean,2)
wls = argrelextrema(cube_2_mean,np.less)
wls = np.asarray(wls)
wls = wls[0]
wls = np.append(0,wls)
print wls
#wls[0] += 15
#wls[3] += 15
N_x_panels = 4
N_y_panels = len(wls)

#x = np.linspace(0,99.0,100)
#x *= 20.8 / 1E3
#y = x
x = np.linspace(0,NX-1,NX)
y = np.linspace(0,NY-1,NY)
x*=20.8*3.0/1E3
y*=20.8*3.0/1E3

#make the size of the figure:
x_size = 2.5
ratio = 0.85
y_size = x_size * float(NY)/float(NX) * ratio

shrinkage = 0.7

noise = 7E-4 * np.sqrt(cube_1_mean[0]*cube_1_mean[:])
residual = (cube1-cube2)
residual[:] /= noise;
residual = residual * residual
print residual.shape
if (int(ifmask)):
	residual[:,:,:,:] *= mask
residual = np.sum(residual,axis=3)
residual = residual[:,:,0] + 0*residual[:,:,3] + 0*residual[:,:,1] + residual[:,:,2]
print 'chisq_max = ', np.amax(residual)
print 'chisq_mean = ', np.mean(residual)
if (int(ifmask)):
	residual /= 4.0 * np.sum(mask)
else:
	residual /= (4.0*NL)
print 'chisq_reduced_max = ', np.amax(residual)
print 'chisq_reduced_mean = ', np.mean(residual)

irange = [0.8,1.2]
vrange = [-1,1]

#N_y_panels -= 2

#plt.figure(figsize=[N_x_panels*x_size,N_y_panels*y_size])
fig, axes = plt.subplots(nrows=N_y_panels,ncols=N_x_panels,figsize=(N_x_panels*x_size,N_y_panels*y_size))
fig.subplots_adjust(right = 0.85,left=0.05,top=0.95,bottom=0.05)
image_no = 0

for j in range (1,N_y_panels+1):
	
	to_plot = np.copy(cube1[:,:,0,wls[j-1]])
	m = np.mean(to_plot)
	to_plot/=m
	s = np.std(to_plot)
	
	#plt.subplot(N_y_panels,N_x_panels,(j-1)*N_x_panels+1)
	ax = axes.flat[image_no]
	im = ax.imshow(to_plot,origin='lower',vmin = irange[0],vmax= irange[1],cmap=cc.cm['fire'],extent=[0,x[-1],0,y[-1]])
	#plt.colorbar(shrink=shrinkage)
	if (j==N_y_panels):
		ax.set_xlabel('$x\,[\mathrm{Mm}]$')
	ax.set_ylabel('$y\,[\mathrm{Mm}]$')
	if (j==1):
		ax.set_title('Observed $I/I_{qs}$')
	image_no+=1

	to_plot = np.copy(cube2[:,:,0,wls[j-1]])
	to_plot/=m
	ax=axes.flat[image_no]
	#plt.subplot(N_y_panels,N_x_panels,(j-1)*N_x_panels+2)
	im = ax.imshow(to_plot,origin='lower',vmin = irange[0],vmax= irange[1],cmap=cc.cm['fire'],extent=[0,x[-1],0,y[-1]])
	#plt.colorbar(shrink=shrinkage)
	if (j==N_y_panels):
		ax.set_xlabel('$x\,[\mathrm{Mm}]$')
	if (j==1):
		ax.set_title('Fitted $I/I_{qs}$')
	image_no+=1

	if (image_no == 2):
		cb_ax = fig.add_axes([0.88, 0.51, 0.03, 0.44])
		cbar = fig.colorbar(im, cax=cb_ax)


	to_plot1 = np.mean(cube1[:,:,3,wls[j-1]:wls[j-1]+6],axis=2)/m*100.0
	s = np.std(to_plot1)

	#plt.subplot(N_y_panels,N_x_panels,(j-1)*N_x_panels+3)
	ax = axes.flat[image_no]
	im = ax.imshow(to_plot1,origin='lower',vmin = vrange[0],vmax=vrange[1],cmap=cc.cm['coolwarm'],extent=[0,x[-1],0,y[-1]])
	#plt.colorbar(shrink=shrinkage)
	if (j==N_y_panels):
		ax.set_xlabel('$x\,[\mathrm{Mm}]$')
	if (j==1):
		ax.set_title('Observed $V/I_{qs}\,[\%]$')
	image_no+=1

	to_plot2 = np.mean(cube2[:,:,3,wls[j-1]:wls[j-1]+6],axis=2)/m*100.0	
	#plt.subplot(N_y_panels,N_x_panels,(j-1)*N_x_panels+4)
	ax = axes.flat[image_no]
	im = ax.imshow(to_plot2,origin='lower',vmin = vrange[0],vmax=vrange[1],cmap=cc.cm['coolwarm'],extent=[0,x[-1],0,y[-1]])
	#plt.colorbar(shrink=shrinkage)
	if (j==N_y_panels):
		ax.set_xlabel('$x\,[\mathrm{Mm}]$')
	if (j==1):
		ax.set_title('Fitted $V/I_{qs}\,[\%]$')
	image_no +=1

	if (image_no == 4):
		cb_ax = fig.add_axes([0.88, 0.05, 0.03, 0.44])
		cbar = fig.colorbar(im, cax=cb_ax)



	

#fig.tight_layout()
#plt.show()
fig.savefig(out_name,fmt='png',bbxox_inches='tight')
fig.savefig(out_name+'.eps',fmt='eps',bbxox_inches='tight')

plt.clf()
plt.cla()
plt.close('all')
