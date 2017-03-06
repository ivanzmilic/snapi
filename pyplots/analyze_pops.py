import numpy as np
import sys

filename = sys.argv[1]

pops = np.loadtxt(filename, unpack = True)
pops = pops.transpose()
pops = pops.reshape(100,30,5)
pops = pops[:,:,2:5]

mean = np.zeros([30,3])
dev = np.zeros([30,3])

for i in range (0,30):
	for j in range (0,3):
		mean[i,j] = np.mean(pops[:,i,j])
		dev[i,j] = np.std(pops[:,i,j])

print dev/mean