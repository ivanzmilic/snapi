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
lambda_l = min(spectrum[:,0])
lambda_m = max(spectrum[:,0])
#spectrum[:,1] /= max(spectrum[:,1])

mpl.rcParams['font.size'] = 18

#plt.ion()
plt.xlabel('$\lambda\,\mathrm{[\AA]}$')
plt.ylabel('$I\,[\mathrm{erg\,cm^{-2}\,s^{-1}}\,\AA^{-1}]$')
plt.xlim([lambda_l, lambda_m])
plt.ylim(min(spectrum[:,1]) - 0.1 * (max(spectrum[:,1]) - min(spectrum[:,1])), max(spectrum[:,1]) * 1.1)
plt.gca().get_xaxis().get_major_formatter().set_useOffset(False)
plt.tight_layout()
plt.savefig(output_file+'_spectrum.eps', format='eps')

#then load numerical and analytical responses

h, wvl, rn = np.loadtxt(rn_file, unpack = True)
h, wvl, ra = np.loadtxt(ra_file, unpack = True)

N_Parameters = 7

rn = rn.reshape(N_Parameters,n_depths, n_wvl)
ra = ra.reshape(N_Parameters,n_depths, n_wvl)

for n in range(0,n_wvl):
	rn[:,:,n] /= spectrum[n,1]
	ra[:,:,n] /= spectrum[n,1]

rn *= 1E4
ra *= 1E4

wvl = wvl.reshape(N_Parameters,n_depths, n_wvl)
h = h.reshape(N_Parameters,n_depths, n_wvl)


wvl = wvl[0][0]
wvl *= 1E8
wvl -= (wvl[n_wvl-1] + wvl[0]) / 2.0
h = h[0,:,0] / 1E5

suffix = ['temperature','density','vt']

for p in range(0,3):

	v_min = np.amin(rn[p])
	v_max = np.amax(rn[p])

	#print v_min, v_max

	#then plot the stuff
	if (rn_file != ra_file):
		plt.clf()
		plt.cla()
		plt.xlim([wvl[0], wvl[n_wvl-1]])
		plt.ylim([h[0], h[n_depths-1]])
		plt.xlabel('$\lambda\,\mathrm{[\AA]}$')
		plt.ylabel('$h\,\mathrm{[km]}$')
		plt.pcolormesh(wvl, h, rn[p], vmin = v_min, vmax = v_max, rasterized=True)
		plt.colorbar()
		plt.tight_layout()
		plt.savefig(output_file+'_numerical_responses_intensity_'+suffix[p]+'.eps', format='eps')

	#v_min = np.amin(ra[p]*100)
	#v_max = np.amax(ra[p]*100)

	#then plot the stuff
	plt.clf()
	plt.cla()
	plt.xlim([wvl[0], wvl[n_wvl-1]])
	plt.ylim([h[0], h[n_depths-1]])
	plt.xlabel('$\lambda\,\mathrm{[\AA]}$')
	plt.ylabel('$h\,\mathrm{[km]}$')
	plt.pcolormesh(wvl, h, ra[p], vmin = v_min, vmax = v_max, rasterized=True)
	plt.colorbar()
	plt.tight_layout()
	plt.savefig(output_file+'_analytical_responses_intensity_'+suffix[p]+'.eps', format='eps')

	rmax = np.amax(np.abs(rn[p]))
	#then plot the relative differences
	if (rn_file != ra_file):
		print np.amax(np.abs(ra[p]-rn[p]))/rmax
		plt.clf()
		plt.cla()
		plt.xlim([wvl[0], wvl[n_wvl-1]])
		plt.ylim([h[0], h[n_depths-1]])
		plt.xlabel('$\lambda\,\mathrm{[\AA]}$')
		plt.ylabel('$h\,\mathrm{[km]}$')
		rel_diff = np.log10(np.abs((ra[p]-rn[p]))/rmax)
		plt.pcolormesh(wvl, h, rel_diff,vmin = -10.0, vmax = -1.0, rasterized=True)
		plt.colorbar()
		plt.tight_layout()
		plt.savefig(output_file+'_relative_difference_responses_intensity_'+suffix[p]+'.eps', format='eps')

plt.clf()
plt.cla()
plt.xlim([wvl[0], wvl[n_wvl-1]])
plt.ylim([h[0], h[n_depths-1]])
plt.xlabel('$\lambda\,\mathrm{[\AA]}$')
plt.ylabel('$h\,\mathrm{[km]}$')
plt.pcolormesh(wvl, h, ra[4], rasterized=True)
plt.colorbar()
plt.tight_layout()
plt.savefig(output_file+'_analytical_responses_intensity_'+'B'+'.eps', format='eps')