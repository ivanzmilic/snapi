import scipy.ndimage.filters as filters
import scipy.signal as flt 
import numpy as np 
import pyana
import sys

file_in = sys.argv[1]
sigma = float(sys.argv[2])

a = pyana.fzread(file_in)
an = a["data"]

dims = an.shape
NP = dims[0]

an_conv = np.zeros(dims)
#an_conv_2 = np.zeros(dims)

for i in range(0,NP):
	an_conv[i] = filters.gaussian_filter(an[i],sigma)
	#an_conv_2[i] = flt.medfilt(an[i],3)

pyana.fzwrite(sys.argv[2]+'_'+file_in,an_conv,0,'bla')
#pyana.fzwrite(sys.argv[2]+'_alt_'+file_in,an_conv_2,0,'bla')