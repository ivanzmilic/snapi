import matplotlib
matplotlib.use('Agg')
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

l = np.loadtxt("lambda_to_fit.dat",unpack=True)
l = l[1]*1E8


for i in range (0,1):
	for j in range(0,1):
		plt.clf()
		plt.cla()
		plt.figure(figsize=[10.0,5.50])

		plt.subplot(221)
		plt.plot(l,obs[i,j,0]/max(obs[i,j,0]),color='red',label='Observation')
		plt.plot(l,fit[i,j,0]/max(obs[i,j,0]),'--',color='blue',label='Fit')
		plt.xlabel('$\lambda\,[\mathrm{\AA}]$')
		plt.ylabel('$\mathrm{Stokes\,I/I_c}$')
		plt.legend()

		plt.subplot(222)
		plt.plot(l,obs[i,j,3]/max(obs[i,j,0]),color='red',label='Observation')
		plt.plot(l,fit[i,j,3]/max(obs[i,j,0]),'--',color='blue',label='Fit')
		plt.xlabel('$\lambda\,[\mathrm{\AA}]$')
		plt.ylabel('$\mathrm{Stokes\,V/I_c}$')

		plt.subplot(223)
		plt.plot(l,obs[i,j,1]/max(obs[i,j,0]),color='red',label='Observation')
		plt.plot(l,fit[i,j,1]/max(obs[i,j,0]),'--',color='blue',label='Fit')
		plt.xlabel('$\lambda\,[\mathrm{\AA}]$')
		plt.ylabel('$\mathrm{Stokes\,Q/I_c}$')
		plt.legend()

		plt.subplot(224)
		plt.plot(l,obs[i,j,2]/max(obs[i,j,0]),color='red',label='Observation')
		plt.plot(l,fit[i,j,2]/max(obs[i,j,0]),'--',color='blue',label='Fit')
		plt.xlabel('$\lambda\,[\mathrm{\AA}]$')
		plt.ylabel('$\mathrm{Stokes\,U/I_c}$')
		
		
		plt.tight_layout()
		plt.savefig('quick_test'+str(i)+'_'+str(j),fmt='png',bbox_inches='tight')
		#plt.savefig('quick_test'+str(i)+'_'+str(j)+'.eps',fmt='eps',bbox_inches='tight')
		plt.close('all')


