import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import numpy as np 
import pyana
import sys
from scipy.signal import argrelextrema
import scipy.ndimage.filters as flt
from matplotlib_scalebar.scalebar import ScaleBar


stokes_cube_file = sys.argv[1]

temp = pyana.fzread(stokes_cube_file)
stokes_cube = temp["data"]
#stokes_cube = np.transpose(stokes_cube,(1,2,0,3))

#compute mean profile and from it deduce the positions of the lines:
mean = np.mean(stokes_cube[:,:,0,:],axis=(0,1))
mean = flt.gaussian_filter(mean,15)
lines=argrelextrema(mean, np.less)
lines=np.asarray(lines[0])
print lines
offset = [0,0,0,5]
NL = stokes_cube.shape[3]
NX = stokes_cube.shape[1]
NY = stokes_cube.shape[0]
N_y = len(lines)
N_x = 4

size = 1.5

#1105

plt.figure(figsize=[4.5*1, 3.7*2])
plt.clf()
plt.cla()
plt.subplot(211)
plt.imshow(stokes_cube[:,:,0,0]/mean[0],origin='lower',cmap='hot')
plt.colorbar(shrink=0.9)
plt.title('$\mathrm{Stokes}\,I\,\mathrm{Fe\,1.56\mu m}$')
#plt.ylabel('$\mathrm{Stokes}\,I$')
plt.subplot(212)
plt.imshow(stokes_cube[:,:,3,370],origin='lower',cmap='coolwarm',vmin=-0.04,vmax=0.04)
plt.title('$V/I$')
plt.colorbar(shrink=0.9)
#plt.tight_layout()
plt.savefig('I_V',bbox_inches='tight')
quit();


plt.clf()
plt.cla()
plt.figure(figsize=[size * float(NX)/float(NY)*N_x,size*N_y])

for l in range(0,N_y):
	for s in range(0,4):
		plt.subplot(N_y,N_x,l*N_x+s+1)
		plt.imshow(stokes_cube[:,:,s,lines[l]+offset[s]],origin='lower')
		plt.colorbar(shrink=0.75)

plt.tight_layout()
plt.savefig(sys.argv[2],fmt='png',bbox_inches='tight')

#quit();
i1 = 150
j1 = 240

i2 = 65
j2 = 130


#plt.clf()
#plt.cla()
#l = np.linspace(15640.0,15670.0,1501)
#plt.plot(l,stokes_cube[i1,j1,0,:]/mean[0])
#plt.plot(l,stokes_cube[i2,j2,0,:]/mean[0])
#plt.xlabel("$\\lambda\,[\mathrm{\AA}]$")
#plt.xlim([l[0],l[-1]])
#plt.ylabel("Normalized $I$")
#plt.savefig('spectra_'+str(i1)+'_'+str(j1),fmt='png',bbox_inches='tight')

l1 = 210
l2 = 210
to_plot=range(210,320)+range(480,530)+range(750,900)


wvl = np.linspace(5887.7914658,5897.18400651,NL)

mean = np.mean(stokes_cube[:,:,0,:],axis=(0,1))
plt.close('all')
counter =0
for l in to_plot:
	counter +=1
	plt.figure(figsize=[8.5,5.5])
	plt.clf()
	plt.cla()
	plt.subplot(311)
	to_show = stokes_cube[:,:,0,l]
	m = np.mean(to_show)
	s = np.std(to_show)
	plt.imshow(stokes_cube[:,:,0,l],origin='lower',cmap='hot',vmin=m-2*s,vmax=m+2*s)
	plt.title('Intensity')
	plt.subplot(312)
	plt.imshow(np.mean(stokes_cube[:,:,3,l-2:l+2],axis=2),origin='lower',vmin=-0.1,vmax=0.1,cmap='coolwarm')
	scalebar = ScaleBar(38.0*1E3) # 1 pixel = 0.2 meter
	plt.gca().add_artist(scalebar)
	plt.tight_layout()
	plt.title('Circular polarization')
	plt.subplot(313)
	plt.plot(wvl,mean/max(mean),color='navy')
	plt.ylim([0,1.1])
	plt.xlim([wvl[0],wvl[-1]])
	plt.xlabel('$\\lambda\,[\mathrm{\AA}]$')
	plt.ylabel('Mean intensity')
	plt.axvline(x=wvl[l],color='navy')
	plt.tight_layout()
	plt.savefig('img_'+str(100+counter),fmt='png',bbox_inches='tight')
	plt.close('all')