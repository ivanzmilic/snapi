import pyana
import numpy as np 
import matplotlib.pyplot as plt 
import sys

ref_in = sys.argv[1]
fit_in = sys.argv[2]

atmos_ref = np.loadtxt(ref_in,unpack=True,skiprows=1)
temp = pyana.fzread(fit_in)
atmos_fit = temp["data"]

atmos_fit = atmos_fit[0,0]

plt.clf()
plt.cla()
plt.title('Temperature comparison')
plt.subplot(211)
plt.plot(atmos_ref[0],atmos_ref[2],label='Reference')
plt.plot(atmos_fit[0],atmos_fit[2],label='Fit')
#plt.plot(atmos_ref[0],atmos_ref[2]-atmos_fit[2],label='Difference')
plt.xlabel('log $\\tau$')
plt.ylabel('T [K]')
plt.legend(loc=4)
plt.subplot(212)
plt.title('LOS velocity comparison')
plt.plot(atmos_ref[0],atmos_ref[9]/1E5,label='Reference')
plt.plot(atmos_fit[0],atmos_fit[9]/1E5,label='Fit')
#plt.plot(atmos_ref[0],atmos_ref[9]/1E5-atmos_fit[9]/1E5,label='Difference')
plt.xlabel('log $\\tau$')
plt.ylabel('V [km/s]')
plt.legend(loc=4)
plt.tight_layout()
plt.savefig('atmos_comp',fmt='png',bbox_inches='tight')

