import numpy as np 
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import sys

spectrum_to_fit = np.loadtxt(sys.argv[1],unpack = True)
spectrum_by_iters = np.loadtxt(sys.argv[2],unpack = True)
#nodes = np.loadtxt(sys.argv[3],unpack = True,skiprows = 1)
chisq = np.loadtxt(sys.argv[4],unpack = True)

spectrum_to_fit[0] *= 1E8

#spectrum_to_fit[0] -= (spectrum_to_fit[0,-1]+spectrum_to_fit[0,0]) * 0.5

N_wlv = spectrum_to_fit.shape[1]
N_iters = spectrum_by_iters.size / (5 * N_wlv)

spectrum_by_iters = spectrum_by_iters.reshape(5,N_iters,N_wlv)

#nodes = nodes.reshape(2,N_iters,4)

limits = np.zeros([2,4])
limits[0,0] = 0.0
limits[1,0] = np.amax(spectrum_by_iters[1,:,:]) * 1.1

for s in range (1,4):
	limits[0,s] = np.amin(spectrum_by_iters[s+1,:,:]) * 1.1
	limits[1,s] = np.amax(spectrum_by_iters[s+1,:,:]) * 1.1

lambda_min = -2.0
lambda_max = 2.0

lambda_min = spectrum_to_fit[0,0]
lambda_max = spectrum_to_fit[0,-1]

# But now also prepare the nodes and their evolution:


#tau_fine = np.linspace(nodes[0,0,0],nodes[0,0,-1],50)
#print N_wlv, N_iters

for i in range (0,N_iters):
	plt.figure(i+1,figsize=(12,6))
	plt.clf()
	plt.cla()
	plt.subplot(231)
	plt.plot(spectrum_to_fit[0], spectrum_to_fit[1],'o')
	plt.plot(spectrum_to_fit[0], spectrum_by_iters[1,i,:],'red',linewidth=2.0)
	plt.ylim([limits[0,0], limits[1,0]])
	plt.xlim([lambda_min,lambda_max])
	plt.ylabel('$\mathrm{Stokes}\,I$')
	plt.subplot(232)
	plt.plot(spectrum_to_fit[0], spectrum_to_fit[2],'o')
	plt.plot(spectrum_to_fit[0], spectrum_by_iters[2,i,:],'red',linewidth=2.0)
	plt.ylim([limits[0,1], limits[1,1]])
	plt.xlim([lambda_min,lambda_max])
	plt.ylabel('$\mathrm{Stokes}\,Q$')

	#plt.subplot(233)
	#f = interp1d(nodes[0,i,:],nodes[1,i,:], kind='cubic')
	#plt.plot(nodes[0,i,:],nodes[1,i,:],'o',markersize=10)
	#plt.plot(tau_fine,f(tau_fine))
	#plt.xlim([-5.5,1.0])
	#plt.ylim([3000,10000])
	#plt.xlabel('$\\tau$')
	#plt.ylabel('$\mathrm{T\,[K]}$')
	
	plt.subplot(234)
	plt.plot(spectrum_to_fit[0], spectrum_to_fit[3],'o')
	plt.plot(spectrum_to_fit[0], spectrum_by_iters[3,i,:],'red',linewidth=2.0)
	plt.ylim([limits[0,2], limits[1,2]])
	plt.xlim([lambda_min,lambda_max])
	plt.ylabel('$\mathrm{Stokes}\,U$')
	plt.xlabel('$\lambda[\AA]$')
	plt.subplot(235)
	plt.plot(spectrum_to_fit[0], spectrum_to_fit[4],'o')
	plt.plot(spectrum_to_fit[0], spectrum_by_iters[4,i,:],'red',linewidth=2.0)
	plt.ylim([limits[0,3], limits[1,3]])
	plt.xlim([lambda_min,lambda_max])
	plt.ylabel('$\mathrm{Stokes}\,V$')
	plt.xlabel('$\lambda[\AA]$')

	plt.subplot(236)

	plt.plot(chisq[0],np.log10(chisq[1]))
	plt.plot(chisq[0,i],np.log10(chisq[1,i]),'o', markersize = 10)
	plt.xlabel('$\mathrm{Iteration\,no.}$')
	plt.ylabel('$\mathrm{log}_{10} \chi^2$')

	plt.tight_layout()
	plt.savefig('fit_'+str(i+1),fmt='png')
