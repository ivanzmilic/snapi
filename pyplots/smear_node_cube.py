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

#Project B, theta, phi -> B_x, B_y B_z

B_z = an[10:14] * np.cos(an[14])
B_t = an[10:14] * np.sin(an[14])

an_conv = np.copy(an)

for i in range(0,NP):
	an_conv[i] = flt.medfilt(an[i],5)
	if (sigma):
		an_conv[i] = filters.gaussian_filter(an[i],sigma)

for i in range(0,0):
	B_z[i] = flt.medfilt(B_z[i],5)
	B_t[i] = flt.medfilt(B_t[i],5)
	if (sigma):
		B_z[i] = filters.gaussian_filter(B_z[i],sigma)
		B_t[i] = filters.gaussian_filter(B_t[i],sigma)

for i in range(0,0):
	an_conv[10+i] = (B_z[i] * B_z[i] + B_t[i] * B_t[i]) ** 0.5

#an_conv[-1] = 0.1


pyana.fzwrite(file_out,an_conv,0,'bla')
