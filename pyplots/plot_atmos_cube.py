import matplotlib.pyplot as plt 
import numpy as np
import sys 
import pyana
import scipy.interpolate.interp1d as interp1d

input_file = sys.argv[1]

temp = pyana.fzread(input_file)
atmos_cube = temp["data"]
print atmos_cube.shape

atmos_cube[:,:,9,:] /= 1E5 #to km/s
atmos_cube[:,:,9,:] *= -1.0 # reverse sign
atmos_cube[:,:,10:12,:] *= 180.0/np.pi #to deg

depths = [35,53,65]
depth_values = atmos_cube[0,0,0,depths]
parameters=[2,9,7,10,11]
param_names = ['T [K]',' v los [km/s]', 'B los [Gauss]', 'theta [deg]', 'phi[deg]']
cmaps=['hot','coolwarm','summer']

#plot a mean stokes profile:
panelsx = len(depths)
panelsy = 3

plt.figure(figsize=[4*panelsx, 3*panelsy])

atmos_cube[:,:,7,:] *= np.cos(atmos_cube[:,:,9,:]*np.pi/180.0)

for p in range(0,panelsy):
	for d in range(0,panelsx):

		
		v_min = np.amin(atmos_cube[:,:,parameters[p],depths[d]])
		v_max = np.amax(atmos_cube[:,:,parameters[p],depths[d]])
		if (p==0 or p==2):
			m = np.mean(atmos_cube[:,:,parameters[p],depths[d]])
			s = np.std(atmos_cube[:,:,parameters[p],depths[d]])
			v_min = m-3*s
			v_max = m+3*s
			if (p==2):
				v_min = 0
				v_max = m+5*s

		plt.subplot(panelsy,panelsx,p*panelsx+d+1)
		plt.imshow(atmos_cube[:,:,parameters[p],depths[d]].transpose(),origin='lower',cmap=cmaps[p],
			vmin=v_min, vmax=v_max)
		plt.colorbar(shrink=0.85)
		if (d==0):
			plt.ylabel(param_names[p])
		if (p==0):
			plt.title('$\log \\tau = $'+str(depth_values[d]))

plt.savefig(input_file+'_atmos.eps',fmt='eps',bbox_inches='tight')
plt.savefig(input_file+'_atmos.png',fmt='png',bbox_inches='tight')