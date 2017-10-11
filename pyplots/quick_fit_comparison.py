import pyana 
import numpy as np 
import matplotlib.pyplot as plt 
import sys

obs_in = sys.argv[1]
fit_in = sys.argv[2]

temp = pyana.fzread(obs_in)
obs = temp["data"]
temp = pyana.fzread(fit_in)
fit = temp["data"]

plt.clf()
plt.cla()
plt.plot(obs[0,0,0],color='red',label='Observation')
plt.plot(fit[0,0,0],color='blue',label='Fit')
plt.legend()
plt.savefig('quick_test',fmt='png')

