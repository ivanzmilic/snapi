import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import numpy as np 
import pyana
import sys
from scipy.signal import argrelextrema
import scipy.ndimage.filters as flt


stokes_cube_file = sys.argv[1]

temp = pyana.fzread(stokes_cube_file)
stokes_cube = temp["data"]

#compute mean profile and from it deduce the positions of the lines:
mean = np.mean(stokes_cube[:,:,0,:],axis=(0,1))
mean = flt.gaussian_filter(mean,5)
lines=argrelextrema(mean, np.less)
lines=np.asarray(lines[0])
print lines
offset = [0,0,0,5]

N_y = len(lines)
N_x = 4

plt.figure(figsize=[12,12])

for l in range(0,N_y):
	for s in range(0,4):
		plt.subplot(N_y,N_x,l*N_x+s+1)
		plt.imshow(stokes_cube[:,:,s,lines[l]+offset[s]],origin='lower')
		plt.colorbar(shrink=0.75)

plt.tight_layout()
plt.savefig('MURAM_ir_lines',fmt='png',bbox_inches='tight')