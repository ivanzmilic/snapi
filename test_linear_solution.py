import numpy as np
import sys

matrix = sys.argv[1]
rhs = sys.argv[2]
point_to_test = int(sys.argv[3])

a = np.loadtxt(matrix, unpack = True)
a = a.transpose()

matrix_size = a.shape[0]

beta = np.loadtxt(rhs, unpack = True)
beta = beta.reshape(-1,matrix_size)

# Do the svd:
U, s, V = np.linalg.svd(a, full_matrices=True)



S = np.zeros([matrix_size, matrix_size])
S = np.diag(s)
# Check if they agree:
print np.allclose(a, np.dot(U, np.dot(S, V)))

#Compute the solution:
s_inv = 1.0 / s
S_inv = np.diag(s_inv)
a_inv = np.dot(V.transpose(), np.dot(S_inv, U.transpose()))
#
print beta.shape
b = beta[point_to_test]
x = np.dot(a_inv, beta[point_to_test])
print x
