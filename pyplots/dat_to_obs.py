import numpy as np 
import pyana
import sys

file_in = sys.argv[1]+'.dat'

spectrum = np.loadtxt(file_in,unpack=True)

NS = 4
NL = spectrum.shape[1]

print 'Number of wavelengths = ', NL

cube = np.zeros([1,1,NS,NL])

cube[0,0] = np.copy(spectrum[1:5])
#for s in range(1,4):
#	cube[0,0,s] /= cube[0,0,0]

pyana.fzwrite(sys.argv[1]+'.f0',cube,0,'placeholder')
