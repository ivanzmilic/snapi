import matplotlib
matplotlib.use('Agg')
import pyana 
import numpy as np 
import matplotlib.pyplot as plt 
import sys

obs_in = sys.argv[1]
fit_in = sys.argv[2]
atmos_in = sys.argv[3]
or_atmos_in = sys.argv[4]
nodes_in = sys.argv[5]

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
nodes = np.transpose(nodes,(0,2,1))

atmos[:,:,7,:] *= np.cos(atmos[:,:,10,:])

model = pyana.fzread(or_atmos_in)["data"]
model[:,:,7,:] *= np.cos(model[:,:,10,:])
model = np.transpose(model,(1,0,2,3))

#fit = np.transpose(fit,(1,0,2,3))

#l = np.linspace(5887.7914658,5897.18400651,967)
#l = l[0:957]
#l = np.linspace(6301,6303,201)
#l = np.linspace(6300.8921,6303.2671,112)
l = np.linspace(15640.0,15670,501)
#l = np.linspace(8540.0,8543.0,151)
#l = l[:151]
l = l[50:451]

x = int(sys.argv[6])
y = int(sys.argv[7])

for i in range (x,x+1):
	for j in range(y,y+1):
		plt.clf()
		plt.cla()
		plt.figure(figsize=[9.0,11.0])

		plt.subplot(421)
		plt.plot(l,obs[i,j,0]/max(obs[i,j,0]),color='red',label='Observation')
		plt.plot(l,fit[i,j,0]/max(obs[i,j,0]),'--',color='blue',label='Fit')
		#plt.plot((fit[i,j,0]/max(obs[i,j,0]) - obs[i,j,0]/max(obs[i,j,0])) * 10.0,label='Residual')
		#plt.xlim([0,22])
		plt.xlabel('$\lambda\,[\mathrm{\AA}]$')
		plt.ylabel('$\mathrm{Stokes\,I/I_c}$')
		#plt.xlim([8541.0,8543.0])
		plt.ylim([0.5,1.2])
		plt.legend()

		plt.subplot(422)
		plt.plot(l,obs[i,j,3]/max(obs[i,j,0]),color='red',label='Observation')
		plt.plot(l,fit[i,j,3]/max(obs[i,j,0]),'--',color='blue',label='Fit')
		#plt.xlim([l[0],l[-1]])
		plt.xlabel('$\lambda\,[\mathrm{\AA}]$')
		plt.ylabel('$\mathrm{Stokes\,V/I_c}$')
		#plt.ylim([0.0,1.4])
		#plt.legend()

		#plt.subplot(423)
		#plt.plot(obs[i,j,1]/max(obs[i,j,0]),color='red',label='Observation')
		#plt.plot(fit[i,j,1]/max(obs[i,j,0]),'--',color='blue',label='Fit')
		#plt.xlim([l[0],l[-1]])
		#plt.xlabel('$\lambda\,[\mathrm{\AA}]$')
		#plt.ylabel('$\mathrm{Stokes\,Q/I_c}$')
		#plt.ylim([0.0,1.4])
		#plt.legend()

		#plt.subplot(424)
		#plt.plot(obs[i,j,2]/max(obs[i,j,0]),color='red',label='Observation')
		#plt.plot(fit[i,j,2]/max(obs[i,j,0]),'--',color='blue',label='Fit')
		#plt.xlim([l[0],l[-1]])
		#plt.xlabel('$\lambda\,[\mathrm{\AA}]$')
		#plt.ylabel('$\mathrm{Stokes\,U/I_c}$')

		#plt.ylim([0.0,1.4])
		#plt.legend()

		plt.subplot(423)
		plt.plot(atmos[i,j,0],atmos[i,j,2],label='Inversion')
		plt.plot(model[i,j,0],model[i,j,2],label='Original MHD atmosphere')
		plt.xlim([-3,1])
		plt.xlabel('$\log\\tau$')
		plt.ylabel('$\mathrm{T\,[K]}$')
		plt.legend()

		plt.subplot(424)
		plt.plot(atmos[i,j,0],atmos[i,j,7])
		plt.plot(model[i,j,0],model[i,j,7])
		plt.xlim([-3,1])
		plt.xlabel('$\log\\tau$')
		plt.ylabel('$B_{\mathrm {los}}\,\mathrm{[Gauss]}$')


		plt.subplot(425)
		plt.plot(atmos[i,j,0],np.log10(atmos[i,j,3]*1.38E-16*atmos[i,j,2]))
		plt.plot(model[i,j,0],np.log10(model[i,j,3]*1.38E-16*model[i,j,2]))
		plt.xlim([-3,1])
		plt.xlabel('$\log\\tau$')
		plt.ylabel('$\log p_{\mathrm {gas}}$')


		plt.subplot(426)
		plt.plot(atmos[i,j,0],-atmos[i,j,9]/1E5)
		plt.plot(model[i,j,0],-model[i,j,9]/1E5)
		plt.xlim([-3,1])
		plt.ylim([-6,6])
		plt.xlabel('$\log\\tau$')
		plt.ylabel('$\mathrm{v_{los}\,[km/s]}$')

		print nodes[:,x,y]



		plt.tight_layout()
		plt.savefig('quick_test'+str(i)+'_'+str(j),fmt='png',bbox_inches='tight')
		plt.savefig('quick_test'+str(i)+'_'+str(j)+'.eps',fmt='eps',bbox_inches='tight')
		plt.close('all')


