import numpy as np 
import matplotlib.pyplot as plt 
import sys

file_in = sys.argv[1]
ND = int(sys.argv[2])
NL = int(sys.argv[3])
file_out = sys.argv[4]

data = np.loadtxt(file_in,unpack=True)

data=data.reshape(5,NL,ND)

h = data[0,0,:] / 1E5 #to km
tau = data[1]
S = data[4]/data[3]

CR = S*np.exp(-tau)*tau

plt.clf()
plt.cla()
plt.plot(h,CR[41]/max(CR[41]))
plt.plot(h,CR[0]/max(CR[0]))
plt.xlim([-100,500])
plt.savefig(file_out,fmt='png')
