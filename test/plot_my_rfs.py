import numpy as np 
import matplotlib.pyplot as plt 
import sys 
#import pyana
# increas the fonts a bit
import matplotlib
matplotlib.rcParams.update({'font.size': 16})

# Load an atmosphere to get the T and N
atmosfile = str(sys.argv[2])
atmos = np.loadtxt(atmosfile, unpack=True, skiprows=1) # might want to load this differently if you have a 3D model
#atmos = pyana.fzread(atmosfile)["data"]
#atmos = atmos[:,0,0,:] # just take the first column

T = atmos[2]
p = atmos[3]
N = p /T / 1.38E-16

filein = sys.argv[1]
NZ = atmos.shape[1]

rf = np.loadtxt(filein, unpack=True)

rf = rf.reshape(6,7,NZ,-1)
print ("info:: the shape of the RF is : ", rf.shape)

# plot the rfs

h = np.copy(rf[0,0,:,0]) / 1E5

ll = np.copy(rf[1,0,0,:]) * 1E8

# The thing is that the RF to temperature is with constant N, but we want with constant pressure, so recall: 
# dI / dT |_P = dI / dT | _N + dI / dN | _T * dN / dT | _P
# P = NkT
# and dN / dT | _P = - N / T

# So to calculate this we need the OG atmosphere:

rf_T_P = rf[2,0,:,:] + rf[2,1,:,:] * (-N[:,None] / T[:,None])

# Then what else we have is we have rf to the particle density, which we do not want, we want to have to the pressure, so:
# dI / dp |_T = dI / dN |_T * dN / dp |_T + dI / dT |_P * dT / dp |_T
# and dN / dp |_T = 1 / kT
# and dT / dp |_T = 0 # is this right? Double check later ;-)
rf_p_T = rf[2,1,:,:] * (1.0 / (1.38E-16 * T[:,None]))

# How to scale the response functions? My favorite is to normalize w.r.t the continuum intensity:
spectrum = np.loadtxt(sys.argv[3], unpack=True)
Icont = spectrum[1,0]

'''plt.figure(figsize=[14,9])
plt.subplot(221)
plt.imshow(rf_T_P/Icont, cmap='PuOr', extent=[ll[0],ll[-1],h[-1],h[0]], aspect='auto',vmin=-3E-5, vmax=3E-5)
plt.colorbar()
plt.title("RF to the temperature with const P")
plt.subplot(222)
plt.imshow(rf[2,0,:,:]/Icont, cmap='PuOr', extent=[ll[0],ll[-1],h[-1],h[0]], aspect='auto',vmin=-3E-5, vmax=3E-5)
plt.colorbar()
plt.title("RF to the temperature with const N")
plt.subplot(223)
plt.imshow(rf[2,3,:,:]/Icont, cmap='bwr', extent=[ll[0],ll[-1],h[-1],h[0]], aspect='auto')
plt.colorbar()
plt.title("RF to the los velocity")
plt.subplot(224)
plt.imshow(rf_p_T/Icont*p[:,None], cmap='inferno', extent=[ll[0],ll[-1],h[-1],h[0]], aspect='auto')
plt.colorbar()
plt.title("RF to the pressure with const T")
plt.tight_layout()
plt.savefig("myrf.png",bbox_inches='tight')'''

# For review paper from SAJ, for NaD lines:
# only 1 x 2 plot, top temperature , below los velocity:
plt.figure(figsize=[10,10])
plt.subplot(211)
plt.imshow(rf[2,0,:,:]/Icont, cmap='plasma', extent=[ll[0],ll[-1],h[-1],h[0]], aspect='auto',vmin=0, vmax=5E-5)
plt.title("RF to the temperature")
plt.xlim(5892.0, 5898.0)
plt.ylim(h[-1],1100)
plt.colorbar()
plt.subplot(212)
plt.imshow(rf[2,3,:,:]/Icont, cmap='bwr', extent=[ll[0],ll[-1],h[-1],h[0]], aspect='auto', vmin = -2E-7, vmax=2E-7)
plt.colorbar()
plt.title("RF to the los velocity")
plt.xlim(5892.0, 5898.0)
plt.ylim(h[-1],1100)
plt.tight_layout()
plt.savefig("myrf.png",bbox_inches='tight')
plt.savefig("Na_D1_rfs_falc.eps",bbox_inches='tight')

# Alsoo plot the spectrum just in case:
plt.clf()
plt.figure(figsize=[8,5])
plt.plot(spectrum[0,:]*1E8,spectrum[1,:])
plt.savefig("spectrum.png",bbox_inches='tight')
