from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import colorcet as cc

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyana
import numpy as np 
import sys
from scipy.interpolate import interp1d
from matplotlib_scalebar.scalebar import ScaleBar


def plot_scatter_hist(x,y):
    xedges, yedges = np.linspace(np.min(x),np.max(x), 42),np.linspace(np.min(x),np.max(y), 42)
    hist, xedges, yedges = np.histogram2d(x, y, (xedges, yedges))
    xidx = np.clip(np.digitize(x, xedges), 0, hist.shape[0]-1)
    yidx = np.clip(np.digitize(y, yedges), 0, hist.shape[1]-1)
    c = hist[xidx, yidx]
    return x,y,c

cube1_in = sys.argv[1]
cube2_in = sys.argv[2]
plot_here = sys.argv[3]

temp = pyana.fzread(cube1_in)
cube1 = temp["data"]
#cube1 = np.transpose(cube1,(1,0,2,3))

temp = pyana.fzread(cube2_in)
cube2 = temp["data"]

print cube1.shape
print cube2.shape
dims = cube1.shape
NX = int(dims[0])
NY = int(dims[1])

x_scale = np.linspace(0,NX-1,NX)
y_scale = np.linspace(0,NY-1,NY)
x_scale *= 20.8/1E3 * 3.0
y_scale *= 20.8/1E3 * 3.0

plt.clf()
plt.cla()
x_panel_size = 3.2*0.8
y_panel_size = 2.0*0.8

c_map = ['','','','']

y_low = 0
y_high = 100
x_low = 0
x_high = 100

tau = [-1.0,0.0,0.5]
tau = np.asarray(tau)
N_tau = tau.size

params = [2,7,9]
maps = ['hot','coolwarm','coolwarm']
suffix = ['T','B','V']
title = ['Temperature [K]', 'LOS B [G]', 'LOS velocity [km/s]']
units = [500.0,100.0,1.0]


cube1[:,:,7,:] *= np.cos(cube1[:,:,10,:])
cube2[:,:,7,:] *= np.cos(cube2[:,:,10,:])
cube1[:,:,9,:] /= -1E5
#cube1[:,:,9,:] -= 0.35
cube2[:,:,9,:] /= -1E5

cube1_to_show = np.zeros([N_tau,NX,NY])
cube2_to_show = np.zeros([N_tau,NX,NY])

scale = [1000.0,1000.0,1.0]

for ii in range (0,3):

	plt.clf()
	plt.cla()
	fig,axes = plt.subplots(N_tau,3,figsize=[3.0*x_panel_size,1.5*len(tau)*y_panel_size])
	fig_no = 0

	p = params[ii]
	
	for i in range(0,NX):
		for j in range(0,NY):
			f = interp1d(cube1[i,j,0,:],cube1[i,j,p,:])
			cube1_to_show[:,i,j] = f(tau)
			f = interp1d(cube2[i+x_low,j+y_low,0,:],cube2[i+x_low,j+y_low,p,:])
			cube2_to_show[:,i,j] = f(tau)


	for i in range (0,N_tau):

		m = np.mean(cube2_to_show[i])
		if (ii>0):
			m = 0
		s = np.std(cube2_to_show[i])
		
		ax = axes.flat[fig_no]
		ax.imshow(cube2_to_show[i],origin='lower',cmap=maps[ii],vmin=m-3*s,vmax=m+3*s,extent=[x_scale[0],x_scale[-1],y_scale[0],y_scale[-1]])
		if (i==0):
			ax.set_title('Simulation')
		ax.set_ylabel("$\log\\tau=$"+"$"+str(tau[i])+"$")
		if (i==N_tau-1):
			ax.set_xlabel('$x\,[\mathrm{Mm}]$')
		print np.mean(cube2_to_show[i])
		fig_no +=1

		m = np.mean(cube1_to_show[i])
		if (ii>0):
			m = 0
		s = np.std(cube1_to_show[i])

		ax = axes.flat[fig_no]
		ax.imshow(cube1_to_show[i],origin='lower',cmap=maps[ii],vmin=m-3*s,vmax=m+3*s,extent=[x_scale[0],x_scale[-1],y_scale[0],y_scale[-1]])
		if (i==0):
			ax.set_title('Inversion')
		print np.mean(cube1_to_show[i])
		if (i==N_tau-1):
			ax.set_xlabel('$x\,[\mathrm{Mm}]$')
		fig_no +=1
		
		ax = axes.flat[fig_no]
		#x,y,c = plot_scatter_hist(cube2_to_show[i].reshape(NX*NY),cube1_to_show[i].reshape(NX*NY))
		#ax.scatter(x,y,c=np.log(c),s=10)
		#ax.plot(cube2_to_show[i].reshape(NX*NY),cube2_to_show[i].reshape(NX*NY),color='red')
		if (i==0):
			ax.set_title(title[ii])
		if (i==N_tau-1):
			ax.set_xlabel('Simulation')
		ax.set_ylabel('Inversion')
		#ax.set_xlim(min(x)-scale[ii],max(x)+scale[ii])
		#ax.set_ylim(min(y)-scale[ii],max(y)+scale[ii])
		fig_no +=1

	for axi in axes.flat:
		axi.xaxis.set_major_locator(plt.MaxNLocator(5))
    	axi.yaxis.set_major_locator(plt.MaxNLocator(5))

	fig.tight_layout()
	fig.savefig(plot_here+'_'+suffix[ii]+'.eps',fmt='eps',bbox_inches='tight')
	fig.savefig(plot_here+'_'+suffix[ii],fmt='png',bbox_inches='tight')
	plt.close('all')