import scipy.ndimage.filters as filters
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


for i in range(0,NP):
	an_conv[i] = filters.gaussian_filter(an[i],sigma)

pyana.fzwrite(sys.argv[2]+'_'+file_in,an_conv,0,'bla')