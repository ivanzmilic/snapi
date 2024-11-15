from matplotlib import rc
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
from scipy.stats.stats import pearsonr


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
transpose = int(sys.argv[4])

temp = pyana.fzread(cube1_in)
cube1 = temp["data"]
if (transpose):
	cube1 = np.transpose(cube1,(1,0,2,3))

temp = pyana.fzread(cube2_in)
cube2 = temp["data"]

print cube1.shape
print cube2.shape
dims = cube1.shape
NX = int(dims[0])
NY = int(dims[1])
dims2 = cube2.shape
NX2 = int(dims2[0])
NY2 = int(dims2[1])

shrink = int(sys.argv[6])
if (shrink):
	NX = NX2
	NY = NY2
x_scale = np.linspace(0,NX2-1,NX2)
y_scale = np.linspace(0,NY2-1,NY2)
x_scale *= 20.8/1E3
y_scale *= 20.8/1E3
scaling  = float(sys.argv[5])


cube1 = np.transpose(cube1,(1,0,2,3))
cube2 = np.transpose(cube2,(1,0,2,3))

plt.clf()
plt.cla()
x_panel_size = 3.2*0.8
y_panel_size = 2.0*0.8

c_map = ['','','','']

tau = [-1.5,-0.5,0.5]
tau = np.asarray(tau)
N_tau = tau.size

params = [2,7,9,8]
maps = ['hot','coolwarm','coolwarm','PuOr']
suffix = ['T','B','V']
title = ['Temperature [K]', 'LOS B [G]', 'LOS velocity [km/s]','turbulent velocity']
units = [500.0,100.0,1.0,1.0]


cube1[:,:,7,:] *= np.cos(cube1[:,:,10,:]) *np.sqrt(4.0*3.141)
cube2[:,:,7,:] *= np.cos(cube2[:,:,10,:]) *np.sqrt(4.0*3.141)
cube1[:,:,9,:] /= -1E5
cube2[:,:,9,:] /= -1E5
cube1[:,:,8,:] /= -1E5
cube2[:,:,8,:] /= -1E5


cube1_to_show = np.zeros([N_tau,NX,NY])
cube2_to_show = np.zeros([N_tau,NX2,NY2])

scale = [1000.0,1000.0,1.0]

for ii in range(0,3):

	plt.clf()
	plt.cla()
	fig,axes = plt.subplots(N_tau,3,figsize=[3.0*x_panel_size,1.5*len(tau)*y_panel_size])
	fig.subplots_adjust(right = 0.85,left=0.05,top=0.95,bottom=0.05)
	fig_no = 0

	print ii
	p = params[ii]
	
	for i in range(0,NX):
		for j in range(0,NY):
			f = interp1d(cube1[i,j,0,:],cube1[i,j,p,:])
			cube1_to_show[:,i,j] = f(tau)
	for i in range(0,NX2):
		for j in range(0,NY2):
			f = interp1d(cube2[i,j,0,:],cube2[i,j,p,:])
			cube2_to_show[:,i,j] = f(tau)


	for i in range (0,N_tau):

		m = np.mean(cube1_to_show[i])
		if (ii>0):
			m = 0
		s = np.std(cube1_to_show[i])
		
		ax = axes.flat[fig_no]
		ax.imshow(cube1_to_show[i],origin='lower',cmap=maps[ii],vmin=m-3*s,vmax=m+3*s,extent=[x_scale[0],x_scale[-1],y_scale[0],y_scale[-1]])
		if (i==0):
			ax.set_title('Original cube')
		ax.set_ylabel("$\log\\tau=$"+"$"+str(tau[i])+"$")
		if (i==N_tau-1):
			ax.set_xlabel('$x\,[\mathrm{Mm}]$')
		else:
			ax.set_xticklabels([])
		#print np.mean(cube2_to_show[i])
		fig_no +=1

		ax = axes.flat[fig_no]
		im=ax.imshow(cube2_to_show[i],origin='lower',cmap=maps[ii],vmin=m-3*s,vmax=m+3*s,extent=[x_scale[0],x_scale[-1],y_scale[0],y_scale[-1]])
		if (i==0):
			ax.set_title('Inversion')
		ax.set_yticklabels([])
		#print np.mean(cube1_to_show[i])
		if (i==N_tau-1):
			ax.set_xlabel('$x\,[\mathrm{Mm}]$')
		else:
			ax.set_xticklabels([])
		
		if (i==0):
			cb_ax = fig.add_axes([0.88, 0.675, 0.03, 0.26])
			#cb_ax = fig.add_axes([0.88, 0.53, 0.03, 0.39])
			cbar = fig.colorbar(im, cax=cb_ax)
		if (i==1):
			cb_ax = fig.add_axes([0.88, 0.37, 0.03, 0.26])
			#cb_ax = fig.add_axes([0.88, 0.075, 0.03, 0.39])
			cbar = fig.colorbar(im, cax=cb_ax)
		if (i==2):
			cb_ax = fig.add_axes([0.88, 0.065, 0.03, 0.26])
			cbar = fig.colorbar(im, cax=cb_ax)
		fig_no +=1
		
		if (NX2 == NX and NY2 == NY):
			ax = axes.flat[fig_no]
		
			ax.hist2d(cube1_to_show[i].flatten(),cube2_to_show[i].flatten(),bins=(100,100),cmap='Purples',vmax=30,range=[[m-3*s,m+3*s],[m-3*s,m+3*s]])
			#print pearsonr(x,y)
			std_difference = np.std(cube1_to_show[i] - cube2_to_show[i])
			#print 'Std of difference = ', std_difference
			mean_mag = np.mean(np.abs(cube2_to_show[i]))
			#print 'Mean magnitude    = ', mean_mag
			#print 'Ratio             = ', std_difference/mean_mag
			ax.plot(cube1_to_show[i].flatten(),cube1_to_show[i].flatten(),color='red')
			if (i==0):
				ax.set_title("Inversion vs Simulation")
			ax.set_xticklabels([])
			ax.set_yticklabels([])
			#if (i==N_tau-1):
			#	ax.set_xlabel('Simulation')
			#ax.set_ylabel('Inversion')
		
			fig_no +=1

	#for axi in axes.flat:
	#	axi.xaxis.set_major_locator(plt.MaxNLocator(6))
    	#axi.yaxis.set_major_locator(plt.MaxNLocator(6))

	#fig.tight_layout()
	fig.subplots_adjust(hspace=0.05, wspace=0.05)
	fig.savefig(plot_here+'_'+suffix[ii]+'.eps',fmt='eps',bbox_inches='tight')
	fig.savefig(plot_here+'_'+suffix[ii],fmt='png',bbox_inches='tight')
	plt.close('all')
