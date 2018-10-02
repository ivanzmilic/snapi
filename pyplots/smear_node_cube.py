import scipy.ndimage.filters as filters
import scipy.signal as flt 
import numpy as np 
import pyana
import sys
import matplotlib.pyplot as plt

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

an_conv = np.copy(an)
start = int(sys.argv[4])
end = int(sys.argv[5])

B_start = start
B_end = end

B_v = an[start:end+1]*np.cos(an[-1])
B_h = an[start:end+1]*np.sin(an[-1])


for i in range(0,B_start):
	#an_conv[i] = flt.medfilt(an[i],5)
	if (sigma):
		an_conv[i] = filters.gaussian_filter(an_conv[i],sigma)

for i in range(B_start,B_end+1):
	if (sigma):
		B_v[i-B_start] = filters.gaussian_filter(B_v[i-B_start],sigma)
		B_h[i-B_start] = filters.gaussian_filter(B_h[i-B_start],sigma)
	an_conv[i] = np.sqrt(B_v[i-B_start] ** 2.0 + B_h[i-B_start] ** 2.0)

an_conv[-1] = np.cos(an_conv[-1])
an_conv[-1] = filters.gaussian_filter(an_conv[-1],sigma)

small = np.where(an_conv[-1] < -0.99)
an_conv[-1,small] = -0.99
big = np.where(an_conv[-1] > 0.99)
an_conv[-1,big] = 0.99
an_conv[-1] = np.arccos(an_conv[-1])

for i in range(0,0):
	plt.clf()
	plt.cla()
	plt.imshow(an_conv[i])
	plt.colorbar()
	plt.savefig('node_'+str(i)+'.png',fmt='png')

pyana.fzwrite(file_out,an_conv,0,'bla')
