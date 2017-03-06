import numpy as np 
import matplotlib.pyplot as plt 
import sys

num = np.loadtxt(sys.argv[1], unpack = True)
ana = np.loadtxt(sys.argv[2], unpack = True)

line1, = plt.plot(ana[0], ana[1])
line2, = plt.plot(num[0], num[1])
plt.legend([line1, line2], ['Analytical', 'Finite Differences'])
plt.ylabel("Log chisq")
plt.xlabel("# of iterations")
plt.yscale('log')
plt.savefig(sys.argv[3], fmt = 'eps')

