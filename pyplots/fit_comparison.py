import matplotlib
matplotlib.use('Agg')
import pyana 
import numpy as np 
import matplotlib.pyplot as plt 
import sys

obs_in = sys.argv[1]
fit1_in = sys.argv[2]
atmos1_in = sys.argv[3]
fit2_in = sys.argv[4]
atmos2_in = sys.argv[5]


temp = pyana.fzread(obs_in)
obs = temp["data"]

temp = pyana.fzread(fit1_in)
fit1 = temp["data"]
temp = pyana.fzread(atmos1_in)
atmos1 = temp["data"]
temp = pyana.fzread(fit2_in)
fit2 = temp["data"]
temp = pyana.fzread(atmos2_in)
atmos2 = temp["data"]

temp = pyana.fzread("inverted_nodes_an.f0")
nodes1 = temp["data"]
temp = pyana.fzread("inverted_nodes_nu.f0")
nodes2 = temp["data"]


temp = pyana.fzread("/home/milic/inversion_code/results/BIFROST_cube/synth/101_300_101_300_Na_full/BIFROST_en024048_bin_1x1_101_300_101_300_atmos.f0")
atmos_mhd = temp["data"]
i_offset = 40
j_offset = 170

l = np.linspace(8540.0,8543.0,151)
l = l[:151]

x = int(sys.argv[6])
y = int(sys.argv[7])

atmos_mhd[:,:,9,:] /= -1E5
atmos1[:,:,9,:] /= -1E5
atmos2[:,:,9,:] /= -1E5

atmos_mhd[:,:,7,:] *= np.cos(atmos_mhd[:,:,10,:])
atmos1[:,:,7,:] *= np.cos(atmos1[:,:,10,:])
atmos2[:,:,7,:] *= np.cos(atmos2[:,:,10,:])

nodes1[13:15,:,:] *= np.cos(nodes1[15,:,:])
nodes2[13:15,:,:] *= np.cos(nodes2[15,:,:])

ana = np.loadtxt("chi_sq_an.dat",unpack=True)
num = np.loadtxt("chi_sq_nu.dat",unpack=True)

err_an = np.loadtxt("errors_an.dat",unpack=True)
err_nu = np.loadtxt("errors_nu.dat",unpack=True)

tau_T = [-6.5,-5.5,-4.5,-3.0,-2.0,-1.0,0.0]
tau_v = [-6.5,-5.0,-3.5,-2.0,-0.5]
tau_B = [-6.0,-0.5]

for i in range (x,x+1):
	for j in range(y,y+1):
		plt.clf()
		plt.cla()
		plt.figure(figsize=[9.0,12.0])

		plt.subplot(421)
		plt.plot(l,obs[i,j,0]/max(obs[i,j,0]),color='red',label='Observation')
		plt.plot(l,fit1[i,j,0]/max(obs[i,j,0]),color='blue',label='Analytical')
		plt.plot(l,fit2[i,j,0]/max(obs[i,j,0]),color='green',label='Finite Differences')
		plt.xlim([l[0],l[-1]])
		plt.xlabel('$\lambda\,[\mathrm{\AA}]$')
		plt.ylabel('$\mathrm{Stokes\,I/I_c}$')
		plt.xlim([8541.0,8543.0])
		plt.ylim([0.0,1.4])
		plt.legend()

		plt.subplot(422)
		plt.plot(l,obs[i,j,3]/max(obs[i,j,0]),color='red',label='Observation')
		plt.plot(l,fit1[i,j,3]/max(obs[i,j,0]),color='blue',label='Analytical')
		plt.plot(l,fit2[i,j,3]/max(obs[i,j,0]),color='green',label='Finite Differences')
		plt.xlim([l[0],l[-1]])
		plt.xlabel('$\lambda\,[\mathrm{\AA}]$')
		plt.ylabel('$\mathrm{Stokes\,V/I_c}$')
		plt.xlim([8541.0,8543.0])
		#plt.legend()

		plt.subplot(423)
		plt.plot(l,obs[i,j,1]/max(obs[i,j,0]),color='red',label='Observation')
		plt.plot(l,fit1[i,j,1]/max(obs[i,j,0]),color='blue',label='Analytical')
		plt.plot(l,fit2[i,j,1]/max(obs[i,j,0]),color='green',label='Finite Differences')
		plt.xlim([l[0],l[-1]])
		plt.xlabel('$\lambda\,[\mathrm{\AA}]$')
		plt.ylabel('$\mathrm{Stokes\,Q/I_c}$')
		plt.xlim([8541.0,8543.0])

		plt.subplot(424)
		plt.plot(l,obs[i,j,2]/max(obs[i,j,0]),color='red',label='Observation')
		plt.plot(l,fit1[i,j,2]/max(obs[i,j,0]),color='blue',label='Analytical')
		plt.plot(l,fit2[i,j,2]/max(obs[i,j,0]),color='green',label='Finite Differences')
		plt.xlim([l[0],l[-1]])
		plt.xlabel('$\lambda\,[\mathrm{\AA}]$')
		plt.ylabel('$\mathrm{Stokes\,U/I_c}$')
		plt.xlim([8541.0,8543.0])

		print err_an[1,0:7]

		plt.subplot(425)
		plt.plot(atmos_mhd[i+i_offset,j+j_offset,0],atmos_mhd[i+i_offset,j+j_offset,2],color='red',label='Original')
		plt.plot(atmos1[i,j,0],atmos1[i,j,2],color='blue',label='Analytical')
		plt.errorbar(tau_T,nodes1[0:7,i,j],yerr=3*err_an[1,0:7],fmt='o',color='blue')
		plt.plot(atmos2[i,j,0],atmos2[i,j,2],color='green',label='Finite Differences')
		plt.errorbar(tau_T,nodes2[0:7,i,j],yerr=3*err_nu[1,0:7],fmt='o',color='green')
		plt.xlabel('$\log \\tau$')
		plt.ylabel('$T\,\mathrm{[K]}$')
		plt.xlim([-6.0,1.0])
		plt.ylim([4000.0,10000])
		plt.legend()

		plt.subplot(426)
		plt.plot(atmos_mhd[i+i_offset,j+j_offset,0],atmos_mhd[i+i_offset,j+j_offset,9],color='red',label='Original')
		plt.plot(atmos1[i,j,0],atmos1[i,j,9],color='blue',label='Analytical')
		plt.errorbar(tau_v,-nodes1[8:13,i,j]/1E5,yerr=3*err_an[1,8:13]/1E5,fmt='o',color='blue')
		plt.plot(atmos2[i,j,0],atmos2[i,j,9],color='green',label='Finite Differences')
		plt.errorbar(tau_v,-nodes2[8:13,i,j]/1E5,yerr=3*err_nu[1,8:13]/1E5,fmt='o',color='green')
		plt.xlabel('$\log \\tau$')
		plt.ylabel('$v_{\mathrm{LOS}}\,\mathrm{[km/s]}$')
		plt.xlim([-6.0,1.0])
		
		plt.subplot(427)
		plt.plot(atmos_mhd[i+i_offset,j+j_offset,0],atmos_mhd[i+i_offset,j+j_offset,7],color='red',label='Original')
		plt.plot(atmos1[i,j,0],atmos1[i,j,7],color='blue',label='Analytical')
		plt.errorbar(tau_B,nodes1[13:15,i,j],yerr=3*err_an[1,13:15],fmt='o',color='blue')
		plt.plot(atmos2[i,j,0],atmos2[i,j,7],color='green',label='Finite Differences')
		plt.errorbar(tau_B,nodes2[13:15,i,j],yerr=3*err_nu[1,13:15],fmt='o',color='green')
		plt.xlabel('$\log \\tau$')
		plt.ylabel('$B_{\mathrm{LOS}}\,\mathrm{[Gauss]}$')
		plt.xlim([-7.0,1.0])

		plt.subplot(428)
		plt.plot(ana[0], ana[1],label='Analytical',color='blue')
		plt.plot(num[0], num[1],label='Finite Differences',color='green')
		#plt.plot(num[0],num[0]/num[0],label="$\chi^2_{r}=1$")
		plt.legend()
		#plt.xlim([1,40])
		plt.ylabel("$\chi^2_{\mathrm{r}}$")
		plt.xlabel("Iteration number")
		plt.yscale('log')
		

		plt.tight_layout()


		#print nodes[:,j,i]

		plt.tight_layout()
		plt.savefig('fit_comparison'+str(i)+'_'+str(j),fmt='png',bbox_inches='tight')
		plt.savefig('fit_comparison'+str(i)+'_'+str(j)+'.eps',fmt='eps',bbox_inches='tight')
		plt.close('all')


