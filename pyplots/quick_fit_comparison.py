import pyana 
import numpy as np 
import matplotlib.pyplot as plt 
import sys

obs_in = sys.argv[1]
fit_in = sys.argv[2]
atmos_in = sys.argv[3]
nodes_in = sys.argv[4]

temp = pyana.fzread(obs_in)
obs = temp["data"]
temp = pyana.fzread(fit_in)
fit = temp["data"]
temp = pyana.fzread(atmos_in)
atmos = temp["data"]
temp = pyana.fzread(nodes_in)
nodes = temp["data"]

l = np.linspace(5887.7914658,5897.18400651,967)
l = l[0:957]

for i in range (136,137):
	for j in range(66,67):
		plt.clf()
		plt.cla()
		plt.figure(figsize=[8.0, 4.5])

		plt.subplot(211)
		plt.plot(l,obs[0,0,0]/max(obs[0,0,0]),color='red',label='Observation')
		plt.plot(l,fit[0,0,0]/max(obs[0,0,0]),color='blue',label='Fit')
		plt.xlim([l[0],l[-1]])
		plt.xlabel('$\lambda\,[\mathrm{\AA}]$')
		plt.ylabel('$\mathrm{Stokes\,I/I_c}$')
		#plt.ylim
		plt.legend()

		plt.subplot(212)
		plt.plot(l,obs[0,0,3]/max(obs[0,0,0]),color='red',label='Observation')
		plt.plot(l,fit[0,0,3]/max(obs[0,0,0]),color='blue',label='Fit')
		plt.xlim([l[0],l[-1]])
		plt.xlabel('$\lambda\,[\mathrm{\AA}]$')
		plt.ylabel('$\mathrm{Stokes\,V/I_c}$')

		#plt.subplot(223)
		#plt.plot(atmos[i,j,0],atmos[i,j,2])
		#plt.xlim([-5,1])
		#plt.xlabel('$\mathrm{log\,}\\tau$')
		#plt.ylabel('$\mathrm{T\,[K]}$')

		#plt.subplot(224)
		#plt.plot(atmos[i,j,0],atmos[i,j,9])
		#plt.xlim([-5,1])
		#plt.xlabel('$\mathrm{log\,}\\tau$')
		#plt.ylabel('$\mathrm{B_{los}\,[Gauss]}$')

		#print nodes[:,j,i]


		plt.tight_layout()
		plt.savefig('quick_test'+str(i)+'_'+str(j),fmt='png',bbox_inches='tight')
		plt.close('all')


