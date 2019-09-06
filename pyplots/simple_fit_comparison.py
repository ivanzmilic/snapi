import pyana 
import matplotlib.pyplot as plt 
import numpy as np 
import sys

obs = pyana.fzread(sys.argv[1])["data"]
fit = pyana.fzread(sys.argv[2])["data"]

print pyana.fzread(sys.argv[3])["data"][:,0,0]

plt.clf()
plt.cla()
plt.subplot(221)
plt.plot(obs[0,0,0])
plt.plot(fit[0,0,0])
plt.subplot(222)
plt.plot(obs[0,0,3])
plt.plot(fit[0,0,3])
plt.subplot(223)
plt.plot(obs[0,0,2])
plt.plot(fit[0,0,2])
plt.subplot(224)
plt.plot(obs[0,0,1])
plt.plot(fit[0,0,1])
plt.tight_layout()
plt.savefig("quick_fit_comparison.png",fmt='png',bbox_inches='tight')