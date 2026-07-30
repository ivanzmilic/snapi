import numpy as np 
import sys

filename = sys.argv[1]

atmosphere = np.loadtxt(filename,unpack=True,skiprows=1)

N_params = atmosphere.shape[0]
N_points = atmosphere.shape[1]


if (N_params == 11):
	N_params += 1
	atmosphere_out = np.zeros([N_params,N_points])
	atmosphere_out[0:N_params-1,:] = atmosphere[:,:]
elif (N_params == 10):
	N_params += 2
	#print atmosphere[8,:]
	atmosphere_out = np.zeros([N_params,N_points])
	atmosphere_out[0:N_params-2,:] = atmosphere[:,:]
else:
	atmosphere_out = np.zeros([N_params,N_points])
	atmosphere_out = atmosphere

#print atmosphere_out[8,:]
#B
atmosphere_out[7,:] = 100. * np.sqrt(2.0)

atmosphere_out[10,:] = np.pi/4.0
#azimuth
atmosphere_out[11,:] = 0.0

atmosphere_out = atmosphere_out.transpose()
np.savetxt(filename, atmosphere_out, fmt = "%2.5e", header = str(N_points)+' '+filename)	