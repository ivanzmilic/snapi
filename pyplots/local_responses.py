import numpy as np 
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
import scipy.ndimage as ndimage

N_parameters = 7

# load two arrays with numerical and analytical population responses
input_n = np.loadtxt(sys.argv[1], dtype = 'double', unpack = True)
input_a = np.loadtxt(sys.argv[2], dtype = 'double', unpack = True)
input_pops = np.loadtxt(sys.argv[3], unpack = True)

output_name = sys.argv[5]

# load populations themselves, but so far we use only the height
atmosphere = np.loadtxt(sys.argv[4], unpack = True, skiprows = 1)

#first two are useless
N_levels = len(input_n[:,0]) - 2
N_depths = int(np.sqrt(len(input_n[0,:])/N_parameters))

input_n  = input_n.reshape(N_levels+2,N_parameters, N_depths, N_depths)
input_a = input_a.reshape(N_levels+2,N_parameters, N_depths, N_depths)

rel_response_n = np.zeros([N_parameters,N_levels, N_depths])
rel_response_a = np.zeros([N_parameters,N_levels, N_depths])

for p in range (0,N_parameters):
	for i in range (2, N_levels+2):
		rel_response_n[p][i-2] = np.diagonal(input_n[i][p]) / input_pops[i]
		rel_response_a[p][i-2] = np.diagonal(input_a[i][p]) / input_pops[i]

h = atmosphere[1] / 1E5

colors = ['Red','Blue','Green','Orange', 'Purple', 'Magenta']
suffix = ['temperature','density', 'vt']

for p in range (0,3):
	for i in range (0, N_levels):
		if (i > 5):
			color_to_plot = colors[5]
		else:
			color_to_plot = colors[i]
		plt.plot(h, rel_response_n[p][i], 'o', color = color_to_plot)
		plt.plot(h, rel_response_a[p][i], '-', color = color_to_plot, label = 'Level '+str(i))

	plt.ylabel('Relative level response')
	plt.xlabel('$h\,[\mathrm{km}]$')
	#plt.ylim([min(rel_response_n[p][i]),max(rel_response_n[p][i])])
	plt.tight_layout()
	plt.legend()
	plt.savefig(output_name+'_local_responses_'+suffix[p]+'.eps', fmt = 'eps')
	plt.clf()
	plt.cla()

	for i in range (0, N_levels-1):
		if (i > 5):
			color_to_plot = colors[5]
		else:
			color_to_plot = colors[i]
		
		plt.plot(h, (rel_response_a[p][i] - rel_response_n[p][i]) / rel_response_n[p][i],'o-', color = color_to_plot, label = 'Level '+str(i))
		
	plt.ylabel('Relative level response - difference')
	plt.xlabel('$h\,[\mathrm{km}]$')
	#plt.ylim([-1,1])
	plt.tight_layout()
	plt.legend()
	plt.savefig(output_name+'_local_responses_difference_'+suffix[p]+'.eps', fmt = 'eps')
	plt.clf()
	plt.cla()


	#now we want to plot populations as a 2d array
	level_to_analyze = int(sys.argv[6])
	v_max = np.ndarray.max(np.abs(input_n[level_to_analyze+2][p] / input_pops[level_to_analyze+2]))
	v_min = np.ndarray.min(input_n[level_to_analyze+2][p] / input_pops[level_to_analyze+2])

	plt.pcolormesh(input_n[level_to_analyze+2][p] / input_pops[level_to_analyze+2], rasterized = True, vmax = v_max, vmin= v_min)
	plt.colorbar()
	plt.tight_layout()
	plt.savefig(output_name+'_level_responses_numerical_'+suffix[p]+'.eps', fmt = 'eps') 
	plt.clf()
	plt.cla()

	plt.pcolormesh(input_a[level_to_analyze+2][p] / input_pops[level_to_analyze+2], rasterized = True, vmax = v_max, vmin= v_min)
	plt.colorbar()
	plt.tight_layout()
	plt.savefig(output_name+'_level_responses_analytical_'+suffix[p]+'.eps', fmt = 'eps') 
	plt.clf()
	plt.cla()

	rel_difference = (input_n[level_to_analyze+2][p] / input_pops[level_to_analyze+2] - input_a[level_to_analyze+2][p] / input_pops[level_to_analyze+2])

	plt.pcolormesh(rel_difference, rasterized = True)
	plt.colorbar()
	plt.tight_layout()
	plt.savefig(output_name+'_level_responses_differences_'+suffix[p]+'.eps', fmt = 'eps') 
	plt.clf()
	plt.cla()

	