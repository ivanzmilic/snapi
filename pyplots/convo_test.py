import numpy as np 
import matplotlib.pyplot as plt 
import sys

def convo_gauss(y,x,sigma,where):
	gauss = 0.56418958354/sigma * np.exp(-(((x[where]-x)/sigma)**2.0))
	to_int = gauss*y
	res = np.trapz(to_int,x)
	return res


array_in = np.loadtxt(sys.argv[1],unpack=True)

sigma = 50.0 #mA
sigma /= 1E11;

N = array_in[0].size
array_mirror = np.zeros(3*N)
array_mirror[0:N] = array_in[1][0]
array_mirror[N:2*N] = array_in[1]
array_mirror[2*N:3*N] = array_in[1][-1]

x_mirror = np.zeros(3*N)
x_mirror[0:N] = array_in[0] - (array_in[0][N-1]-array_in[0][0])
x_mirror[N:2*N] = array_in[0]
x_mirror[2*N:3*N] = array_in[0] + (array_in[0][N-1]-array_in[0][0])

plt.clf()
plt.plot(x_mirror,array_mirror)
plt.savefig("iwantutoshowme.png",fmt='png')


result = np.zeros(N)

for i in range(0,N):
	result[i] = convo_gauss(array_mirror,x_mirror,sigma,i+N)

plt.clf()
plt.plot(array_in[0]*1E8,array_in[1])
plt.plot(array_in[0]*1E8,result)
plt.savefig("convo_testing.png",fmt='png')




