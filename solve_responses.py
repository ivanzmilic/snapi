import numpy as np 
import sys

matrix_file = sys.argv[1]
rhs_file = sys.argv[2]

matrix = np.loadtxt(matrix_file,unpack = True)
matrix = matrix.transpose()
rhs = np.loadtxt(rhs_file, unpack = True)

ND = int(sys.argv[3])
NL = int(sys.argv[4])
rhs = rhs.reshape(ND,ND*NL)

U,w,V = np.linalg.svd(matrix, full_matrices = True)

pinv = np.linalg.inv(matrix)
test = np.dot(pinv,matrix)
#print test.diagonal()
np.savetxt("test1.txt", test)

w_inv = np.diag(1.0/w)

pinv = np.dot(np.dot(V.T,w_inv),U.T)

np.savetxt("test2.txt", np.dot(pinv,matrix))


response = np.zeros([ND,ND*NL])

for l in range(0,ND):
	response[l] = np.dot(pinv,rhs[l])

response = response.reshape(ND*ND,NL)

response_wo = np.zeros([7*ND*ND,NL+2])
response_wo[0:ND*ND,2:NL+2] = response[:,:]

np.savetxt("responses_python.txt", response_wo)


