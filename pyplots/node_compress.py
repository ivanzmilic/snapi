import pyana
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import numpy as np
import pywt
import sys
import skimage

from skimage.restoration import (denoise_wavelet, estimate_sigma)
from skimage.measure import compare_psnr


nodes_in = sys.argv[1] #nodes file
percentage = float(sys.argv[2]) # compression factor

#read in the nodes
temp = pyana.fzread(nodes_in)
nodes = temp["data"]
NP = nodes.shape[0]
print 'Total number of parameters = ',NP
NX = nodes.shape[1]
NY = nodes.shape[2]

#however here we try some other form of wavelet denoising:
#norm = np.amax(np.abs(nodes[index]))
#nodes[index]/=norm
#sigma_est = estimate_sigma(nodes[index], multichannel=True, average_sigmas=True)
#print sigma_est
#reconstructed = denoise_wavelet(nodes[index],wavelet='db8', multichannel=True,method='VisuShrink', mode='hard',sigma=5.0*sigma_est)
#reconstructed *= norm
#nodes[index] *=norm
nodes[NP-2] *= np.cos(nodes[NP-1])
nodes[NP-1,:,:] = 0.1
for p in range(0,NP-1):

	#do 2D wavelet transformation
	coeffs = pywt.dwt2(nodes[p],'db8')

	#assign it to these mighty arrays
	cA,(cH,cV,cD) = coeffs
	#print cA.shape
	#regularize
	cH[:,:] = 0.0
	cV[:,:] = 0.0
	cD[:,:] = 0.0
	#cA[47:57,:] = 0.0

	#reconstruct
	reconstructed = pywt.idwt2(coeffs,'db8')

	plt.clf()
	plt.cla()
	plt.subplot(311)
	plt.imshow(nodes[p])
	plt.colorbar()
	plt.subplot(312)
	plt.imshow(reconstructed)
	plt.colorbar()
	plt.subplot(313)
	plt.imshow(nodes[p]-reconstructed)
	plt.colorbar()
	plt.tight_layout()
	plt.savefig('test'+str(p),fmt='png',bbox_inches='tight')
	plt.close('all')

	nodes[p] = reconstructed

pyana.fzwrite(sys.argv[3],nodes,0,'placeholder')



