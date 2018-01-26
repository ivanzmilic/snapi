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

#print fit[0,0,0]
#print obs[0,0,0]
temp = pyana.fzread(atmos_in)
atmos = temp["data"]
atmos.shape
temp = pyana.fzread(nodes_in)
nodes = temp["data"]

fit = np.transpose(fit,(1,0,2,3))

#l = np.linspace(5887.7914658,5897.18400651,967)
#l = l[0:957]
#l = np.linspace(6301,6303,201)
#l = np.linspace(6300.8921,6303.2671,112)

x = int(sys.argv[5])
y = int(sys.argv[6])

for i in range (x,x+1):
	for j in range(y,y+1):
		plt.clf()
		plt.cla()
		plt.figure(figsize=[5.5,18.0])

		plt.subplot(611)
		plt.plot(obs[i,j,0]/max(obs[i,j,0]),color='red',label='Observation')
		plt.plot(fit[i,j,0]/max(obs[i,j,0]),color='blue',label='Fit')
		#plt.xlim([l[0],l[-1]])
		plt.xlabel('$\lambda\,[\mathrm{\AA}]$')
		plt.ylabel('$\mathrm{Stokes\,I/I_c}$')
		plt.ylim([0.0,1.4])
		plt.legend()

		plt.subplot(612)
		plt.plot(obs[i,j,3]/max(obs[i,j,0]),color='red',label='Observation')
		plt.plot(fit[i,j,3]/max(obs[i,j,0]),color='blue',label='Fit')
		#plt.xlim([l[0],l[-1]])
		plt.xlabel('$\lambda\,[\mathrm{\AA}]$')
		plt.ylabel('$\mathrm{Stokes\,V/I_c}$')

		plt.subplot(613)
		plt.plot(atmos[i,j,0],atmos[i,j,2])
		plt.xlim([-5,1])
		plt.xlabel('$\mathrm{log\,}\\tau$')
		plt.ylabel('$\mathrm{T\,[K]}$')

		plt.subplot(614)
		plt.plot(atmos[i,j,0],atmos[i,j,7]*np.cos(atmos[i,j,10]))
		plt.xlim([-5,1])
		plt.xlabel('$\mathrm{log\,}\\tau$')
		plt.ylabel('$\mathrm{B_{los}\,[Gauss]}$')
		
		plt.subplot(615)
		plt.plot(atmos[i,j,0],atmos[i,j,9]/1E5)
		plt.xlim([-5,1])
		plt.xlabel('$\mathrm{log\,}\\tau$')
		plt.ylabel('$v_{los}$')

		plt.subplot(616)
		plt.plot(atmos[i,j,0],atmos[i,j,8]/1E5)
		plt.xlim([-5,1])
		plt.xlabel('$\mathrm{log\,}\\tau$')
		plt.ylabel('$v_{turb}$')

		print nodes[:,j,i]

		plt.tight_layout()
		plt.savefig('quick_test'+str(i)+'_'+str(j),fmt='png',bbox_inches='tight')
		plt.close('all')


