import matplotlib.pyplot as plt 
import numpy as np
import sys 
import pyana
from matplotlib_scalebar.scalebar import ScaleBar

input_file = sys.argv[1]

temp = pyana.fzread(input_file)
atmos_cube = temp["data"]
print atmos_cube.shape

atmos_cube[:,:,9,:] /= 1E5 #to km/s
atmos_cube[:,:,9,:] *= -1.0 # reverse sign
#atmos_cube[:,:,10:12,:] *= 180.0/np.pi #to deg


depths = [16,26,34]
depth_values = atmos_cube[0,0,0,depths]
parameters=[2,9]
param_names = ['T [K]','B los [Gauss]','v los [km/s]','theta [deg]', 'phi[deg]']
cmaps=['hot','coolwarm','coolwarm']

#plot mean stokes profile:
panelsx = len(depths)
panelsy = len(parameters)

print depth_values

#plt.clf()
#plt.cla()
#i = 220
#plt.figure(figsize=[6.5,3.5])
#plt.imshow(atmos_cube[i,:,2,:].transpose(),vmin=4000.0,vmax=10000.0,cmap='coolwarm',extent=[0,288*20.6/1E3,-800,600],aspect='auto')
#plt.xlabel('$x\,\mathrm{[Mm]}$')
#plt.ylabel('$z\,\mathrm{[Km]}$')
#plt.title('Vertical temperature structure')
#plt.colorbar()

#plt.savefig('vertical_slice.png',fmt='png',bbox_inches='tight')
#plt.close('all')
plt.clf()
plt.cla()
#quit();


plt.figure(figsize=[4.5*panelsx, 3.7*panelsy])
#plt.figure(figsize=[8,3.5])

atmos_cube[:,:,7,:] *= np.cos(atmos_cube[:,:,10,:])

for p in range(0,panelsy):
	for d in range(0,panelsx):

		
		v_min = np.amin(atmos_cube[:,:,parameters[p],depths[d]])
		v_max = np.amax(atmos_cube[:,:,parameters[p],depths[d]])
		if (p==0 or p==1):
			m = np.mean(atmos_cube[:,:,parameters[p],depths[d]])
			s = np.std(atmos_cube[:,:,parameters[p],depths[d]])
			if (p>0):
				m=0
			v_min = m-3*s
			v_max = m+3*s
			#if (p==2):
			#	v_min = 0
			#	v_max = m+5*s

		plt.subplot(panelsy,panelsx,p*panelsx+d+1)
		plt.imshow(atmos_cube[:,:,parameters[p],depths[d]].transpose(),origin='lower',cmap=cmaps[p],
			vmin=v_min, vmax=v_max)
		plt.colorbar(shrink=0.7)
		if (d==0):
			plt.title(param_names[p])
		#if (p==0):
			#plt.title('$\mathrm{Photosphere}$')

		#if (p>0):
			#scalebar = ScaleBar(20.8*1E3) # 1 pixel = 0.2 meter
			#plt.gca().add_artist(scalebar)
		#plt.tight_layout()

plt.savefig(input_file+'_atmos.eps',fmt='eps',bbox_inches='tight')
plt.savefig(input_file+'_atmos.png',fmt='png',bbox_inches='tight')