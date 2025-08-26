import numpy as np 
import matplotlib.pyplot as plt 
import sys 

filein = sys.argv[1]
NZ = int(sys.argv[2])
NL = int(sys.argv[3])

rf = np.loadtxt(filein, unpack=True)

rf = rf.reshape(6,7,NZ,NL)
print ("info:: the shape of the RF is : ", rf.shape)

# plot the rfs

h = np.copy(rf[0,0,:,0]) / 1E5

ll = np.copy(rf[1,0,0,:]) * 1E8

plt.figure(figsize=[14,9])
plt.subplot(221)
plt.imshow(rf[2,0,:,:], cmap='inferno', extent=[ll[0],ll[-1],h[-1],h[0]], aspect='auto')
plt.colorbar()
plt.title("RF to the temperature")
plt.subplot(222)
plt.imshow(rf[2,1,:,:], cmap='cividis', extent=[ll[0],ll[-1],h[-1],h[0]], aspect='auto')
plt.colorbar()
plt.title("RF to the particle density")
plt.subplot(223)
plt.imshow(rf[2,3,:,:], cmap='bwr', extent=[ll[0],ll[-1],h[-1],h[0]], aspect='auto')
plt.colorbar()
plt.title("RF to the los velocity")
plt.subplot(224)
plt.imshow(rf[2,1,:,:], cmap='inferno', extent=[ll[0],ll[-1],h[-1],h[0]], aspect='auto')
plt.colorbar()
plt.title("RF to the pressure - not yet there - ignore")
plt.tight_layout()
plt.savefig("myrf.png",bbox_inches='tight')