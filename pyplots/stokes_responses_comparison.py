import numpy as np 
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt 
import scipy.ndimage as ndimage

#turn input arguments into something usable
profile_file = sys.argv[1]
rn_file = sys.argv[2]
ra_file = sys.argv[3]
n_depths = int(sys.argv[4])
n_wvl = int(sys.argv[5])
output_file = sys.argv[6]

#first load and plot the profile

spectrum = np.loadtxt(profile_file)
plt.plot(spectrum[:,0] * 1E8, spectrum[:,1])
spectrum[:,0] *= 1E8
spectrum[:,0] -= (spectrum[-1,0] + spectrum[0,0]) * 0.5 

lambda_l = min(spectrum[:,0]+0.2)
lambda_m = max(spectrum[:,0]-0.2)
#spectrum[:,1] /= max(spectrum[:,1])

mpl.rcParams['font.size'] = 10


#plt.ion()
plt.figure(1)
plt.clf()
plt.cla()
plt.subplot(221)
plt.xlabel('$\lambda\,\mathrm{[\AA]}$')
plt.ylabel('$I\,[\mathrm{erg\,cm^{-2}\,s^{-1}}\,\AA^{-1}]$')
plt.xlim([lambda_l, lambda_m])
plt.ylim(min(spectrum[:,1]) - 0.1 * (max(spectrum[:,1]) - min(spectrum[:,1])), max(spectrum[:,1]) * 1.1)
plt.plot(spectrum[:,0], spectrum[:,1])
plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
plt.tight_layout()
plt.subplot(222)
plt.xlabel('$\lambda\,\mathrm{[\AA]}$')
plt.ylabel('$Q\,[\mathrm{erg\,cm^{-2}\,s^{-1}}\,\AA^{-1}]$')
plt.xlim([lambda_l, lambda_m])
plt.plot(spectrum[:,0], spectrum[:,2])
plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
plt.tight_layout()
plt.subplot(223)
plt.xlabel('$\lambda\,\mathrm{[\AA]}$')
plt.ylabel('$U\,[\mathrm{erg\,cm^{-2}\,s^{-1}}\,\AA^{-1}]$')
plt.xlim([lambda_l, lambda_m])
plt.plot(spectrum[:,0], spectrum[:,3])
plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
plt.tight_layout()
plt.subplot(224)
plt.xlabel('$\lambda\,\mathrm{[\AA]}$')
plt.ylabel('$V\,[\mathrm{erg\,cm^{-2}\,s^{-1}}\,\AA^{-1}]$')
plt.xlim([lambda_l, lambda_m])
plt.plot(spectrum[:,0], spectrum[:,4])
plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
plt.tight_layout()

plt.savefig(output_file+'_stokes_spectrum.eps', format='eps')

#then load numerical and analytical responses

inputn = np.loadtxt(rn_file, unpack = True)
inputa = np.loadtxt(ra_file, unpack = True)

h = np.copy(inputn[0])
wvl = np.copy(inputn[1])
rn = np.copy(inputn[2:6])
ra = np.copy(inputa[2:6])

N_Parameters = 7
N_Parameters_to_plot = 7

rn = rn.reshape(4,N_Parameters,n_depths, n_wvl)
ra = ra.reshape(4,N_Parameters,n_depths, n_wvl)

for n in range(0,n_wvl):
	rn[:,:,:,n] /= spectrum[n,1]
	ra[:,:,:,n] /= spectrum[n,1]

rn *= 1E4
ra *= 1E4

wvl = wvl.reshape(N_Parameters,n_depths, n_wvl)
h = h.reshape(N_Parameters,n_depths, n_wvl)

wvl = wvl[0][0]
wvl *= 1E8
wvl -= (wvl[n_wvl-1] + wvl[0]) / 2.0

suffix = ['temperature','density','vt','vmacro','B', 'theta', 'phi']

h = h[0,:,0]
#h/= 1E5

hmax = -5
hmin = h[-1]

#yname = '$\log\,\\tau_{500}$'
yname = '$h\,[\mathrm{km}]$'
for p in range(0,7):

	v_min = np.zeros(4)
	v_max = np.zeros(4)
	for s in range(0,4):
		v_min[s] = np.amin(rn[s][p])
		v_max[s] = np.amax(rn[s][p])
		
	#print v_min, v_max

	#then plot the stuff
	if (rn_file != ra_file):
		plt.figure(1)
		plt.clf()
		plt.cla()
		plt.subplot(221)
		plt.xlim([lambda_l, lambda_m])
		plt.ylim([hmin, hmax])
		plt.title('$\mathrm{Stokes}\,I$')
		#plt.xlabel('$\lambda\,\mathrm{[\AA]}$')
		plt.ylabel(yname)
		plt.pcolormesh(wvl, h, rn[0,p,:,:], vmin = v_min[0], vmax = v_max[0], rasterized=True)
		plt.colorbar()
		plt.tight_layout()
		plt.subplot(222)
		plt.xlim([lambda_l, lambda_m])
		plt.ylim([hmin, hmax])
		plt.title('$\mathrm{Stokes}\,Q$')
		#plt.xlabel('$\lambda\,\mathrm{[\AA]}$')
		#plt.ylabel(yname)
		plt.pcolormesh(wvl, h, rn[1,p,:,:], vmin = v_min[1], vmax = v_max[1], rasterized=True)
		plt.colorbar()
		plt.tight_layout()
		plt.subplot(223)
		plt.title('$\mathrm{Stokes}\,U$')
		plt.xlim([lambda_l, lambda_m])
		plt.ylim([hmin, hmax])
		plt.xlabel('$\lambda\,\mathrm{[\AA]}$')
		#plt.ylabel(yname)
		plt.pcolormesh(wvl, h, rn[2,p,:,:], vmin = v_min[2], vmax = v_max[2], rasterized=True)
		plt.colorbar()
		plt.tight_layout()
		plt.subplot(224)
		plt.xlim([lambda_l, lambda_m])
		plt.ylim([hmin, hmax])
		plt.title('$\mathrm{Stokes}\,V$')
		plt.xlabel('$\lambda\,\mathrm{[\AA]}$')
		plt.ylabel(yname)
		plt.pcolormesh(wvl, h, rn[3,p,:,:], vmin = v_min[3], vmax = v_max[3], rasterized=True)
		plt.colorbar()
		plt.tight_layout()
		plt.savefig(output_file+'_numerical_responses_intensity_'+suffix[p]+'.eps', format='eps')

	#v_min = np.amin(ra[p]*100)
	#v_max = np.amax(ra[p]*100)

	#then plot the stuff
	plt.figure(1);
	plt.clf()
	plt.cla()
	plt.subplot(221)
	plt.xlim([lambda_l, lambda_m])
	plt.ylim([hmin, hmax])
	plt.title('$\mathrm{Stokes}\,I$')
	#plt.xlabel('$\lambda\,\mathrm{[\AA]}$')
	plt.ylabel(yname)
	plt.pcolormesh(wvl, h, ra[0,p,:,:], vmin = v_min[0], vmax = v_max[0], rasterized=True)
	plt.colorbar()
	plt.tight_layout()
	plt.subplot(222)
	plt.xlim([lambda_l, lambda_m])
	plt.ylim([hmin, hmax])
	plt.title('$\mathrm{Stokes}\,Q$')
	#plt.xlabel('$\lambda\,\mathrm{[\AA]}$')
	#plt.ylabel(yname)
	plt.pcolormesh(wvl, h, ra[1,p,:,:], vmin = v_min[1], vmax = v_max[1], rasterized=True)
	plt.colorbar()
	plt.tight_layout()
	plt.subplot(223)
	plt.xlim([lambda_l, lambda_m])
	plt.ylim([hmin, hmax])
	plt.title('$\mathrm{Stokes}\,U$')
	plt.xlabel('$\lambda\,\mathrm{[\AA]}$')
	plt.ylabel(yname)
	plt.pcolormesh(wvl, h, ra[2,p,:,:], vmin = v_min[2], vmax = v_max[2], rasterized=True)
	plt.colorbar()
	plt.tight_layout()
	plt.subplot(224)
	plt.xlim([lambda_l, lambda_m])
	plt.ylim([hmin, hmax])
	plt.title('$\mathrm{Stokes}\,V$')
	plt.xlabel('$\lambda\,\mathrm{[\AA]}$')
	#plt.ylabel(yname)
	plt.pcolormesh(wvl, h, ra[3,p,:,:], vmin = v_min[3], vmax = v_max[3], rasterized=True)
	plt.colorbar()
	plt.tight_layout()
	plt.savefig(output_file+'_analytical_responses_intensity_'+suffix[p]+'.eps', format='eps')

	rmax = np.zeros(4)
	for s in range(0,4): rmax[s] = np.amax(np.abs(rn[s][p]))
	#then plot the relative differences, only for stokes I
	if (rn_file != ra_file):
		rel_diff = np.log10(np.abs((ra[0,p,:,:]-rn[0,p,:,:]))/rmax[0])
		#print 10.0**np.amax(rel_diff)
		plt.clf()
		plt.cla()
		plt.xlim([lambda_l, lambda_m])
		plt.ylim([h[0], h[n_depths-1]])
		plt.xlabel('$\lambda\,\mathrm{[\AA]}$')
		plt.ylabel('$h\,\mathrm{[km]}$')
		plt.pcolormesh(wvl, h, rel_diff,vmin = -10.0,vmax=-1.0, rasterized=True)
		plt.colorbar()
		plt.tight_layout()
		plt.savefig(output_file+'_relative_difference_responses_intensity_'+suffix[p]+'.eps', format='eps')
		print suffix[p]
		for s in range(0,4):
			rel_diff = (np.abs((ra[s,p,:,:]-rn[s,p,:,:]))/rmax[s])
			print np.amax(rel_diff)

