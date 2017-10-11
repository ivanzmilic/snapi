import matplotlib.pyplot as plt 
import numpy as np
import sys 
import pyana

input_file = sys.argv[1]

temp = pyana.fzread(input_file)
Stokes_cube = temp["data"]
print Stokes_cube.shape

#plot a mean stokes profile:
I_mean = np.mean(Stokes_cube,axis=(0,1))
#NL = 501
NL = 501
l = np.linspace(5892,5897,NL)
#l = np.linspace(5887,5893.3,201)
l/= 1E8
S = np.zeros([5,NL])
S[0] = l
S[1:5] = I_mean
np.savetxt('BIFROST_Na_mean_profile.dat',S.transpose(),fmt='%1.7e')

plt.clf()
plt.cla()
plt.plot(I_mean[0])
plt.savefig('mean_spectrum.png',fmt='png')

#print the images for the movie:

panelsx = 4
panelsy = 4

plt.figure(figsize=[4*panelsx, 4*panelsy])

wvl = np.array([265,280,287,295])
wvl_c = 295
sname=['Stokes I','Q/I','U/I','V/I']

for s in range(1,4):
	Stokes_cube[:,:,s,:] /= I_mean[0]
	Stokes_cube[:,:,s,:] *= 100.0 #to percent

mydpi = 80

for l in range(251,350):
	plt.clf()
	plt.cla()
	plt.figure(figsize=(366/float(mydpi),1680/float(mydpi)),dpi=mydpi)
	plt.subplot(511)
	to_plot = Stokes_cube[:,:,0,l]
	m = np.mean(to_plot)
	s = np.std(to_plot)
	plt.imshow(to_plot,origin='lower',cmap='hot',vmin=m-3*s,vmax=m+3*s)
	plt.subplot(512)
	to_plot = Stokes_cube[:,:,1,l]
	plt.imshow(to_plot,origin='lower',cmap='coolwarm',vmin=-0.5,vmax=0.5)
	plt.subplot(513)
	to_plot = Stokes_cube[:,:,2,l]
	plt.imshow(to_plot,origin='lower',cmap='coolwarm',vmin=-0.5,vmax=0.5)
	plt.subplot(514)
	to_plot = Stokes_cube[:,:,3,l]
	plt.imshow(to_plot,origin='lower',cmap='coolwarm',vmin=-5,vmax=5)
	plt.subplot(515)
	plt.plot(S[1])
	plt.axvline(x=l,color='red')
	plt.xlim([50,450])
	#plt.tight_layout()

	plt.savefig('i_'+str(l+100),fmt='png')
	plt.close('all')




for s in range(0,4):
	for l in range(0,4):

		to_plot = 0
		v_min = 1
		v_max = 1
		if (s==0):
			m = np.mean(Stokes_cube[:,:,s,wvl[l]])
			to_plot = Stokes_cube[:,:,s,wvl[l]] / m
			m = np.mean(to_plot)
			std = np.std(to_plot)
			
			v_min = m-4*std
			v_max = m+4*std
			if (v_min<0):
				v_min=0
		elif (s==3):
			to_plot = Stokes_cube[:,:,s,wvl[l]]
			v_min = -5
			v_max = 5
		else:
			to_plot = Stokes_cube[:,:,s,wvl[l]]
			v_min = -0.5
			v_max = 0.5
			
		plt.subplot(panelsy,panelsx,s*panelsx+l+1)
		plt.imshow(to_plot,vmin=v_min,vmax=v_max,cmap='gray',origin='lower')
		#
		#if (l==4):
		#	plt.colorbar(shrink=0.5)
		#plt.tight_layout()
		if (s==0):
			plt.title('$'+str((wvl_c-wvl[l])*10)+'\,\mathrm{m \AA}$')
		if (l==0):
			plt.ylabel(sname[s])
		
plt.savefig(input_file+'_plot.eps',fmt='eps',bbox_inches='tight')
plt.savefig(input_file+'_plot.png',fmt='png',bbox_inches='tight')






