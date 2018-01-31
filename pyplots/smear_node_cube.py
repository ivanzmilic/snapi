import scipy.ndimage.filters as filters
import scipy.signal as flt 
import numpy as np 
import pyana
import sys

def svd_compress(A,limit):
	
	U,s,V = np.linalg.svd(A,full_matrices=True)

	dims = A.shape

	small = np.where(s<max(s)*limit)
	s[small] = 0.0
	s_m = np.zeros(dims)
	print small
	s_m[:dims[0],:dims[0]] = np.diag(s)

	A_sparse = np.dot(U,np.dot(s_m,V))

	return A_sparse

file_in = sys.argv[1]
sigma = float(sys.argv[2])
file_out = sys.argv[3]

a = pyana.fzread(file_in)
an = a["data"]

dims = an.shape
NP = dims[0]

an_conv = np.zeros(dims)

#an[-1] = np.cos(an[-1])
#an[8:11] *= an[-1]
#an[-1] = np.arccos(an[-1])
#an[8:11] /= np.cos(an[-1])

#an[8:10] *= np.sum(an[8:10])/np.sum(an[10])

an_conv = np.copy(an)

for i in range(0,NP-5):
	an_conv[i] = flt.medfilt(an[i],7)
	an_conv[i] = filters.gaussian_filter(an[i],sigma)

pyana.fzwrite(file_out,an_conv,0,'bla')
