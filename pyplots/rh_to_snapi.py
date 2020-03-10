#transforms an rh atmosphere to snapi (spinor) atmosphere
# one day I should make Rh-like input to the snapi eh :-) 

import sys
import numpy as np 

file_in = sys.argv[1]
ND = int(sys.argv[2])
file_out = sys.argv[3]

# you need to hard code skips
skip1 = 12
skip2 = 106
in1 = np.loadtxt(file_in,unpack=True,skiprows=skip1,max_rows=ND)
in2 = np.loadtxt(file_in,unpack=True,skiprows=skip2,max_rows=ND)

atmos = np.zeros([12,ND])

#height
atmos[1] = in1[0] * 1E5 #km to cm
atmos[2] = in1[1]
atmos[3] = (np.sum(in2,axis=0)*1.1+in1[2]) * 1.38E-16 * atmos[2]
atmos[7] = 0.0 # usually no mag field
atmos[8] = in1[4]*1E5 #km/s to cm/2
atmos[9] = in1[3]*1E5 # same 

np.savetxt(file_out,atmos.T,fmt="%1.5e")
