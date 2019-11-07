from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=False)

import matplotlib
matplotlib.use('Agg')
import pyana
import numpy as np 
import matplotlib.pyplot as plt 
import sys
from scipy.signal import argrelextrema
import scipy.ndimage.filters as flt
from matplotlib import ticker
from astropy.io import fits

cube1_in = sys.argv[1]
cube2_in = sys.argv[2]
out_name = sys.argv[3]
ifmask = sys.argv[4]
maskfile = sys.argv[5]
fmt = int(sys.argv[6])

if (fmt == 1):
	temp = pyana.fzread(cube1_in)
	cube1 = temp["data"]
	#cube1 = np.transpose(cube1,(1,0,2,3))
	temp = pyana.fzread(cube2_in)
	cube2 = temp["data"]

else:
	cube1 = fits.open(cube1_in)[0].data 
	cube2 = fits.open(cube2_in)[0].data 


if (int(ifmask)):
	mask = np.loadtxt(maskfile,skiprows=1)
	print mask.shape


dims = cube1.shape
NX = dims[0]
NY = dims[1]
NL = dims[3]
print NX, NY

#before plotting the spatial distributions, plot averaged spectra:

cube_1_mean = np.mean(cube1[:,:,0,:],axis=(0,1))
cube_2_mean = np.mean(cube2[:,:,0,:],axis=(0,1))

cube_1_mean = flt.gaussian_filter(cube_1_mean,5)
wls = argrelextrema(cube_1_mean,np.less)
wls = np.asarray(wls)
wls = wls[0]
wls = np.append(0,wls)
print wls
N_x_panels = 4
N_y_panels = len(wls)

plt.plot(cube_1_mean)
plt.savefig('test.png')
plt.clf()
plt.cla()

x = np.linspace(0,NX-1,NX)
y = np.linspace(0,NY-1,NY)
x*=20.8/1E3
y*=20.8/1E3

#make the size of the figure:
x_size = 2.5
ratio = 0.85
y_size = x_size * float(NY)/float(NX) * ratio

shrinkage = 0.7

irange = [0.8,1.2]
vrange = [-3,3]

shift = 2
#restrict ourselves only to single wavelengths:
to_plot_1 = np.zeros([NX,NY,4,len(wls)])
to_plot_1[:,:,0,:] = cube1[:,:,0,wls]
to_plot_1[:,:,3,:] = cube1[:,:,3,wls+shift]
del	cube1
cube1=1.0
to_plot_2 = np.zeros([NX,NY,4,len(wls)])
to_plot_2[:,:,0,:] = cube2[:,:,0,wls]
to_plot_2[:,:,3,:] = cube2[:,:,3,wls+shift]
del cube2
cube2=1.0


fig, axes = plt.subplots(nrows=N_y_panels,ncols=N_x_panels,figsize=(N_x_panels*x_size,N_y_panels*y_size))
fig.subplots_adjust(right = 0.85,left=0.05,top=0.95,bottom=0.1)
image_no = 0



for j in range (1,N_y_panels+1):
	
	to_plot = np.copy(to_plot_1[:,:,0,j-1])
	m = np.mean(to_plot)
	to_plot/=m
	s = np.std(to_plot)
	
	#Observed intensity:
	ax = axes.flat[image_no]
	im = ax.imshow(to_plot,origin='lower',vmin = irange[0],vmax= irange[1],cmap='hot',extent=[0,x[-1],0,y[-1]])
	if (j==N_y_panels):
		ax.set_xlabel('$x\,[\mathrm{Mm}]$')
	else:
		ax.set_xticklabels([])
	ax.set_ylabel('$y\,[\mathrm{Mm}]$')
	if (j==1):
		ax.set_title('Observed $I/I_{qs}$')
	image_no+=1

	#Fitted intensity:
	to_plot = np.copy(to_plot_2[:,:,0,j-1])
	to_plot/=m
	ax=axes.flat[image_no]
	im = ax.imshow(to_plot,origin='lower',vmin = irange[0],vmax= irange[1],cmap='hot',extent=[0,x[-1],0,y[-1]])
	if (j==N_y_panels):
		ax.set_xlabel('$x\,[\mathrm{Mm}]$')
	if (j!=N_y_panels):
		ax.set_xticklabels([])
	if (j==1):
		ax.set_title('Fitted $I/I_{qs}$')
	ax.set_yticklabels([])
	image_no+=1

	if (image_no == 2):
		cb_ax = fig.add_axes([0.88, 0.55, 0.03, 0.4])
		cbar = fig.colorbar(im, cax=cb_ax)

	#Observed V/I_qs
	to_plot1 = to_plot_1[:,:,3,j-1]/m*100.0
	s = np.std(to_plot1)

	ax = axes.flat[image_no]
	im = ax.imshow(to_plot1*np.sqrt(4.0*3.14),origin='lower',vmin = vrange[0],vmax=vrange[1],cmap='coolwarm',extent=[0,x[-1],0,y[-1]])
	if (j==N_y_panels):
		ax.set_xlabel('$x\,[\mathrm{Mm}]$')
	if (j!=N_y_panels):
		ax.set_xticklabels([])
	if (j==1):
		ax.set_title('Observed $V/I_{qs}\,[\%]$')
	ax.set_yticklabels([])
	image_no+=1

	#Fitted V/I_qs
	to_plot2 = to_plot_2[:,:,3,j-1]/m*100.0
	
	ax = axes.flat[image_no]
	im = ax.imshow(to_plot2*np.sqrt(4.0*3.14),origin='lower',vmin = vrange[0],vmax=vrange[1],cmap='coolwarm',extent=[0,x[-1],0,y[-1]])
	if (j==N_y_panels):
		ax.set_xlabel('$x\,[\mathrm{Mm}]$')
	if (j!=N_y_panels):
		ax.set_xticklabels([])
	if (j==1):
		ax.set_title('Fitted $V/I_{qs}\,[\%]$')
	ax.set_yticklabels([])
	image_no +=1

	if (image_no == 4):
		cb_ax = fig.add_axes([0.88, 0.1, 0.03, 0.4])
		cbar = fig.colorbar(im, cax=cb_ax)

fig.subplots_adjust(hspace=0.05, wspace=0.05)
fig.savefig(out_name,fmt='png',bbxox_inches='tight')
fig.savefig(out_name+'.eps',fmt='eps',bbxox_inches='tight')

plt.clf()
plt.cla()
plt.close('all')
