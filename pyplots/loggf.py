import numpy as np

A = 6E5

llambda = 5892.0E-10 # in m

g_l = 1
g_u = 3

c = 2.997E8
me = 9.1E-31
e = 1.6E-19
eps_0 = 8.854E-12

nu = c/llambda

f12 = A * eps_0 * me * c**3.0 * g_u / g_l / e**2.0 / 2.0 / np.pi / nu **2.0

print (np.log10(g_l * f12))