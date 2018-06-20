import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import numpy as np 
import sys
import pyana

file_in = sys.argv[1]
file_out = sys.argv[2]

a = pyana.fzread(file_in)
cube = a["data"]

#first transform Q,U,V to physical units:
for s in range(1,4):
	cube[:,:,s,:] *= cube[:,:,0,:]

spectrum_mean = np.mean(cube,axis=(0,1))

#find median values:
zero_level = np.zeros(4)
for s in range (0,4):
	zero_level[s] = np.median(spectrum_mean[s])

#correct for I-> crosstalk
for s in range (1,4):
	spectrum_mean[s] -= spectrum_mean[0] * zero_level[s]/zero_level[0]

#ad-hoc correction for V->Q,U xtalk

spectrum_mean[3] += spectrum_mean[1]+spectrum_mean[2]
spectrum_mean[1] -= spectrum_mean[1]
spectrum_mean[2] -= spectrum_mean[2]

#but also for the whole cube:
for s in range (1,4):
	cube[:,:,s,:] -= cube[:,:,0,:] * zero_level[s]/zero_level[0]

cube[:,:,3,:] += cube[:,:,1,:]+cube[:,:,2,:]
cube[:,:,1,:] = cube[:,:,2,:] = 0

for s in range (1,4):
	cube[:,:,s,:] /= cube[:,:,0,:]

wherenan = np.isnan(cube)
cube[wherenan] = 0.0

spectrum_to_plot=spectrum_mean
plt.subplot(221)
plt.plot(spectrum_to_plot[0])
plt.subplot(222)
plt.plot(spectrum_to_plot[1])
plt.subplot(223)
plt.plot(spectrum_to_plot[2])
plt.subplot(224)
plt.plot(spectrum_to_plot[3])
plt.tight_layout()
plt.savefig('mean_spectrum.png',fmt='png')

pyana.fzwrite(file_out,cube,0,'bla')
