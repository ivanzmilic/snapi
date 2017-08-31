import numpy as np 
import scipy.signal as sgn
import matplotlib.pyplot as plt 
import sys 
import pyana 
import matplotlib

input_fitted = sys.argv[1] #fitted cube
input_lambdaf = sys.argv[2] #lambda used for fitting

input_obs = sys.argv[3] #observations we have fitted, not necessary same dimension as fitted
input_lambdao = sys.argv[4] #lambda of original observations, same as above

input_nodes = sys.argv[5] #inferred values of the nodes
input_atmos = sys.argv[6] #inferred atmospheres

print_maps_here = sys.argv[7]

#mask = np.loadtxt(sys.argv[7])

#matplotlib.rcParams['figure.figsize'] = 7, 10 #do I even want this?

#l_fit = np.loadtxt(input_lambdaf,unpack = True)
#l_obs = np.loadtxt(input_lambdao,unpack = True)

#l_fit[1]*=1E8
#l_obs[1]*=1E8

a = pyana.fzread(input_fitted)
fitted_cube = a["data"]
dims = fitted_cube.shape

l_offset = 650

#Here we read the parameter map.
pin = pyana.fzread(input_nodes)
parameters = pin["data"]
temp = parameters.shape
NN = temp[0] #total number of nodes

#keep in mind this one is transposed:
NY = dims[1]
NX = dims[0]
NL = dims[3]

b = pyana.fzread(input_obs)
obs_cube = b["data"]

#print obs_cube[-1,:,0,837-l_offset+15]
#print fitted_cube[:,-1,0,837-l_offset+15]
#print parameters[0,:,-1]


a_read = pyana.fzread(input_atmos)
atmospheres = a_read["data"]

#l = np.loadtxt(input_lambdaf,unpack=True)
#wvl = l[1]*1E8

#print atmospheres.shape

#print parameters[:,0,0]
#print parameters[5,0,0]*np.cos(parameters[7,0,0])
#print parameters[6,0,0]*np.cos(parameters[7,0,0])

V_weak_field = np.copy(fitted_cube[:,:,3,:])

for i in range(0,NX):
	for j in range(0,NY):
		V_weak_field[i,j,:] = np.gradient(fitted_cube[i,j,0,:])/0.00972312703583
		B_los = 0.0#parameters[5,i,j]*np.cos(parameters[4,i,j])
		V_weak_field[i,j,:] *= -4.697E-13*1.33*(5896.0**2.0)*B_los

V_total_obs = np.zeros([NX,NY])
V_total_fit = np.zeros([NX,NY])

for i in range(0,NX):
	for j in range (0,NY):
		line_center = np.where(fitted_cube[i,j,0] == min(fitted_cube[i,j,0]))
		lc = line_center[0]
		V_flipped = np.copy(fitted_cube[i,j,3])
		V_flipped[lc:] *= -1.0
		V_total_fit[i,j] = np.sum(V_flipped)
		V_flipped = np.copy(obs_cube[j,i,3])
		V_flipped[lc:] *= -1.0
		V_total_obs[i,j] = np.sum(V_flipped)

#Prepare chisq:
chisq = np.zeros([NX,NY])

noise = 3E12

for i in range(0,NX):
	for j in range (0,NY):
		chisq[i,j] = np.sum(((fitted_cube[i,j,0,:]-obs_cube[j,i,0,:])/noise)**2.0)
		chisq[i,j] += 4.0*np.sum(((fitted_cube[i,j,3,:]-obs_cube[j,i,3,:])/noise)**2.0)
		chisq[i,j] /= (NL-9.0)
		



#Here we pring out some profiles
xl=0
xh=0
yl=0
yh=0
for i in range(xl,xh+1):
	for j in range (yl,yh+1):
		plt.clf()
		plt.cla()
		plt.figure(figsize=[12,10])
		plt.subplot(321)
		plt.plot(fitted_cube[i,j,0,:])
		plt.plot(obs_cube[j,i,0,:])
		plt.xlabel("Wavelength")
		plt.ylabel("Stokes I")
		
		plt.subplot(322)

		plt.plot(fitted_cube[i,j,3,:],label='Fitted')
		plt.plot(obs_cube[j,i,3,:],label='Observed')

		#plt.plot(V_weak_field[i][j],label='WF from obs')
		
		plt.xlabel("Wavelength")
		plt.ylabel("Stokes V")

		plt.subplot(323)
		plt.plot(atmospheres[i,j,0],atmospheres[i,j,1])
		plt.xlabel("$\log \\tau$")
		plt.ylabel("$\mathrm{T\,[K]}$")
		plt.subplot(324)
		plt.plot(atmospheres[i,j,0],atmospheres[i,j,3]/1E5)
		plt.xlabel("$\log \\tau$")
		plt.ylabel("$v_{los}$")
		plt.subplot(325)
		plt.plot(atmospheres[i,j,0],atmospheres[i,j,4])
		plt.xlabel("$\log \\tau$")
		plt.ylabel("B [Gauss]")
		
		plt.tight_layout()
		plt.savefig("test_cube"+str(i)+"_"+str(j)+".png")
		plt.close('all')


		
if (NX<=1 and NY<=1):
	print parameters[:,0,0]
	quit();

#Hard-coded plotting of some images.

barshrink = 0.8

plt.clf()
plt.cla()

#-----------------------------------------------------------------------------------------------------
# NOW NODES THEMSELVES -------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------

T_nodes = [0,1,2,3]
T_nodes_tau = [-2.4,-1.6,-0.8,0.0]
#vt_nodes = [5]
vs_nodes = [4,5,6]
vs_nodes_tau = [-3.5,-2.0,-0.5]
B_nodes = [7,8,9]
B_nodes_tau = [-3.5,-2.0,-0.5]
theta_nodes = [10]

panelsx=4
panelsy=7

#pre-determined wavelengths
l_core_Na = 837-l_offset
l_core_Fe = 505-l_offset
l_core_Ni = 523-l_offset
l_c       = 655-l_offset

Tmap = 'hot'
Vmap = 'coolwarm'
Bmap = 'coolwarm'
Imap = 'hot'
Pmap = 'Spectral'
Dmap = 'coolwarm'

plt.clf()
plt.cla()

defsize = 3
plt.figure(figsize=[defsize*3*panelsx, defsize*panelsy])

intstart = 4 #row where plotting the intensity starts

#Nodes panels:
#Row1&2: Temperature and microturbulentce

for i in range(T_nodes[0],T_nodes[-1]+1):
	plt.subplot(panelsy,panelsx,i+1)
	plt.imshow(parameters[i],origin='lower',cmap=Tmap)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.title('Temperature at $\log\,\\tau$ = '+str(T_nodes_tau[i-T_nodes[0]]))


#for i in range(vt_nodes[0],vt_nodes[-1]+1):
#	parameters[i] /= 1E5
#	plt.subplot(panelsy,panelsx,i+1)
#	plt.imshow(parameters[i],origin='lower')
#	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
#	plt.title('Microturbulent velocity')

for i in range(vs_nodes[0],vs_nodes[-1]+1):

	parameters[i] /= -1E5
	m = np.mean(parameters[i])
	s = np.std(parameters[i])

	plt.subplot(panelsy,panelsx,i+1)
	plt.imshow(parameters[i],origin='lower',cmap=Vmap,vmin=-3*s,vmax=3*s)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.title('Systematic velocity at $\log\,\\tau$ = '+str(vs_nodes_tau[i-vs_nodes[0]]))

s = np.copy(B_nodes)
for i in range(B_nodes[0],B_nodes[-1]+1):
	parameters[i] *= np.cos(parameters[theta_nodes[0]])
	s[i-B_nodes[0]] = 3.0*np.std(parameters[i])


for i in range(B_nodes[0],B_nodes[-1]+1):
	plt.subplot(panelsy,panelsx,i+2)
	plt.imshow(parameters[i],origin='lower',cmap=Bmap,vmin=-1500,vmax=1500)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.title('$\mathrm{B\,[Gauss]}$ at $\log\,\\tau =$'+str(B_nodes_tau[i-B_nodes[0]]))

#m = np.mean(chisq)
#s = np.std(chisq)
#plt.subplot(panelsy,panelsx,B_nodes[-1]+2)
#plt.imshow(chisq,origin='lower',cmap=Tmap,vmin=0,vmax=m+5*s)
#plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
#plt.title('$\chi^2_{\mathrm{reduced}}$')

#Add B gradient if possible:

#B_grad = (parameters[B_nodes[0]]-parameters[B_nodes[1]])
#m = np.mean(B_grad)
#s = np.std(B_grad)
#plt.subplot(panelsy,panelsx,8)
#plt.imshow(B_grad,origin='lower',cmap='coolwarm',vmin=m-3*s,vmax=m+3*s)
#plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
#plt.title('$\mathrm{B\,[Gauss]}$ difference')


# Right hand side panels:
#Continuum intensity
if (l_c >= 0):

	i_conto = np.copy(obs_cube[:,:,0,l_c].transpose())
	i_c_mean = np.mean(i_conto)
	i_conto /= i_c_mean
	sigma = np.std(i_conto)
	
	plt.subplot(panelsy,panelsx,intstart*panelsx+1)
	plt.imshow(i_conto,origin='lower',vmin = 1.0-3*sigma,vmax = 1.0+3*sigma,cmap=Imap)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.title('Observed continuum intensity')

	i_contf = np.copy(fitted_cube[:,:,0,l_c])
	i_contf /= i_c_mean

	plt.subplot(panelsy,panelsx,intstart*panelsx+2)
	plt.imshow(i_contf,origin='lower',vmin = 1.0-3*sigma,vmax = 1.0+3*sigma,cmap=Imap)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.title('Fitted continuum intensity')

	plt.subplot(panelsy,panelsx,(intstart)*panelsx+3)
	plt.imshow((i_conto-i_contf)/i_conto,origin='lower',cmap=Dmap,vmin=-0.2,vmax=0.2)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.title('Continuum intensity difference')


#Fe intensity
if (l_core_Fe >= 0):

	i_conto = np.copy(obs_cube[:,:,0,l_core_Fe].transpose())
	i_c_mean = np.mean(i_conto)
	i_conto /= i_c_mean
	sigma = np.std(i_conto)

	plt.subplot(panelsy,panelsx,intstart*panelsx+3)
	plt.imshow(i_conto,origin='lower',vmin = 1.0-3*sigma,vmax = 1.0+3*sigma,cmap=Imap)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.title('Observed Fe line core')

	i_contf = np.copy(fitted_cube[:,:,0,l_core_Fe])
	i_contf /= i_c_mean

	plt.subplot(panelsy,panelsx,intstart*panelsx+4)
	plt.imshow(i_contf,origin='lower',vmin = 1.0-3*sigma,vmax = 1.0+3*sigma,cmap=Imap)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.title('Fitted Fe line core')

	plt.subplot(panelsy,panelsx,(intstart+1)*panelsx+2)
	plt.imshow((i_conto-i_contf)/i_conto,origin='lower',cmap=Dmap,vmin=-0.2,vmax=0.2)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.title('Fe intensity difference')


#Ni intensity
if (l_core_Ni >= 0):

	i_conto = np.copy(obs_cube[:,:,0,l_core_Ni].transpose())
	i_c_mean = np.mean(i_conto)
	i_conto /= i_c_mean
	sigma = np.std(i_conto)

	plt.subplot(panelsy,panelsx,intstart*panelsx+5)
	plt.imshow(i_conto,origin='lower',vmin = 1.0-3*sigma,vmax = 1.0+3*sigma,cmap=Imap)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.title('Observed Ni line core')

	i_contf = np.copy(fitted_cube[:,:,0,l_core_Ni])
	i_contf /= i_c_mean

	plt.subplot(panelsy,panelsx,intstart*panelsx+6)
	plt.imshow(i_contf,origin='lower',vmin = 1.0-3*sigma,vmax = 1.0+3*sigma,cmap=Imap)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.title('Fitted Ni line core')

	plt.subplot(panelsy,panelsx,(intstart+1)*panelsx+3)
	plt.imshow((i_conto-i_contf)/i_conto,origin='lower',cmap=Dmap,vmin=-0.2,vmax=0.2)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.title('Ni intensity difference')


#Na intensity
if (l_core_Na >= 0):

	i_conto = np.copy(obs_cube[:,:,0,l_core_Na+15].transpose())
	i_c_mean = np.mean(i_conto)
	i_conto /= i_c_mean
	sigma = np.std(i_conto)

	#print i_conto[:,-1]

	plt.subplot(panelsy,panelsx,intstart*panelsx+5)
	plt.imshow(i_conto,origin='lower',vmin = 1.0-3*sigma,vmax = 1.0+3*sigma,cmap=Imap)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.xlim([0,NY-1])
	plt.ylim([0,NX-1])
	plt.title('Observed Na D1 line core')

	i_contf = np.copy(fitted_cube[:,:,0,l_core_Na+15])
	i_contf /= i_c_mean

	plt.subplot(panelsy,panelsx,intstart*panelsx+6)
	plt.imshow(i_contf,origin='lower',vmin = 1.0-3*sigma,vmax = 1.0+3*sigma,cmap=Imap)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.xlim([0,NY-1])
	plt.ylim([0,NX-1])
	plt.title('Fitted Na D1 line core')

	plt.subplot(panelsy,panelsx,intstart*panelsx+7)
	plt.imshow((i_conto-i_contf)/i_conto,origin='lower',cmap=Dmap,vmin=-0.2,vmax=0.2)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.xlim([0,NY-1])
	plt.ylim([0,NX-1])
	plt.title('Na D1 intensity difference')

#STOKES SIGNAL:
V_shift = -15
if (l_core_Na >= 0):

	V = obs_cube[:,:,3,l_core_Na+V_shift].transpose()/obs_cube[:,:,0,l_c].transpose()
	V[np.isnan(V)] = 0
	
	u = 5*np.std(V)
	d = -5*np.std(V)
	
	plt.subplot(panelsy,panelsx,(intstart+2)*panelsx+1)
	plt.imshow(V,origin='lower',cmap=Pmap,vmin=d, vmax=u)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.xlim([0,NY-1])
	plt.ylim([0,NX-1])
	plt.title('Observed Na Stokes V')

	V = fitted_cube[:,:,3,l_core_Na+V_shift]/fitted_cube[:,:,0,l_c]
	
	plt.subplot(panelsy,panelsx,(intstart+2)*panelsx+2)
	plt.imshow(V,origin='lower',cmap=Pmap,vmin=d, vmax=u)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.xlim([0,NY-1])
	plt.ylim([0,NX-1])
	plt.title('Fitted Na Stokes V')

	#V = V_weak_field[:,:,l_core_Na+V_shift]/obs_cube[:,:,0,l_c]
	
	#plt.subplot(panelsy,panelsx,(intstart+1)*panelsx+6)
	#plt.imshow(V,origin='lower',cmap=Pmap,vmin=d, vmax=u)
	#plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	#plt.xlim([0,NY-1])
	#plt.ylim([0,NX-1])
	#plt.title('Weak field Na Stokes V')


	u = 5*np.std(V_total_obs)
	d = -5*np.std(V_total_obs)
	plt.subplot(panelsy,panelsx,(intstart+2)*panelsx+3)
	plt.imshow(V_total_obs,origin='lower',cmap=Pmap,vmin=d, vmax=u)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.xlim([0,NY-1])
	plt.ylim([0,NX-1])
	plt.title('Total Stokes V - obs')

	plt.subplot(panelsy,panelsx,(intstart+2)*panelsx+4)
	plt.imshow(V_total_fit,origin='lower',cmap=Pmap,vmin=d, vmax=u)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.xlim([0,NY-1])
	plt.ylim([0,NX-1])
	plt.title('Total Stokes V - fit')

	rel_diff = (V_total_fit-V_total_obs)/np.amin(V_total_obs)
	print np.std(rel_diff) 
	plt.subplot(panelsy,panelsx,(intstart+1)*panelsx+4)
	plt.imshow(rel_diff,origin='lower',cmap=Pmap,vmin=-0.5, vmax=0.5)
	plt.colorbar(fraction=0.046, pad=0.04,shrink=barshrink)
	plt.xlim([0,NY-1])
	plt.ylim([0,NX-1])
	plt.title('Total Stokes V - relative difference')





plt.tight_layout()
plt.savefig(print_maps_here+'.eps',fmt='eps')
plt.savefig(print_maps_here+'.png',fmt='png')
plt.clf()
plt.cla()















