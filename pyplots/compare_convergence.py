import numpy as np 
import matplotlib.pyplot as plt 
import sys

num = np.loadtxt(sys.argv[1], unpack = True)
ana = np.loadtxt(sys.argv[2], unpack = True)

line1, = plt.plot(ana[0], ana[1])
line2, = plt.plot(num[0], num[1])
line3, = plt.plot(num[0],num[0]/num[0])
plt.legend([line1, line2, line3], ['Analytical', 'Finite Differences',"$\chi^2_{r}=1$"])
plt.xlim([1,40])
plt.ylabel("$\log \chi^2$")
plt.xlabel("Iteration number")
plt.yscale('log')

plt.savefig(sys.argv[3]+'.eps', fmt = 'eps')
plt.savefig(sys.argv[3], fmt = 'png')

