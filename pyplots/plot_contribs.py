import numpy as np 
import matplotlib.pyplot as plt

atmosphere_size = 10

#load numerical and analytical ones
index_pert, index_response, n1n, n2n, n3n = np.loadtxt('../responses_numerical.txt', unpack=True)

index_pert, index_response, n1a, n2a, n3a = np.loadtxt('../responses_analytical.txt', unpack=True)

#load falc, just to have referent height

tau_continuum, h, T = np.loadtxt('../cfg/debug.dat', skiprows = 1, unpack=True, usecols=(0,1,2))

#re-order response
n1n = n1n.reshape(atmosphere_size, atmosphere_size)
n2n = n2n.reshape(atmosphere_size, atmosphere_size)
n3n = n3n.reshape(atmosphere_size, atmosphere_size)
#n4n = n4n.reshape(atmosphere_size, atmosphere_size)
n1a = n1a.reshape(atmosphere_size, atmosphere_size)
n2a = n2a.reshape(atmosphere_size, atmosphere_size)
n3a = n3a.reshape(atmosphere_size, atmosphere_size)
#n4a = n4a.reshape(atmosphere_size, atmosphere_size)

#now start plotting.

h/= 1E5 #go to km

plt.xlabel("Height [km]")
plt.ylabel("Relative response [%]")

#plt.plot(h, n1n.diagonal() * 100, 'b-')
#plt.plot(h, n1a.diagonal() * 100, 'r-')


#plt.plot(h, n2n.diagonal() * 100, 'b-.')
#plt.plot(h, n2a.diagonal() * 100, 'r-.')

plt.plot(h, n3n.diagonal() * 100, 'b--')
plt.plot(h, n3a.diagonal() * 100, 'r--')


#plt.plot(h, n4n.diagonal() * 100, 'b:')
#plt.plot(h, n4a.diagonal() * 100, 'r:')

#plt.axis([1900, 2100, -1, 1])

#plt.show()
#plt.clf()

plt.savefig('num_vs_an.png', format='png')



