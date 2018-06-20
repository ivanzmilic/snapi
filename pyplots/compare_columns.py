import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyana
import numpy as np 
import sys

temp = pyana.fzread(sys.argv[1])
cube_1 = temp["data"]
temp = pyana.fzread(sys.argv[2])
cube_2 = temp["data"]

#offsets, applied to the first cube only
x_low = int(sys.argv[3])
x_high = int(sys.argv[4])
y_low = int(sys.argv[5])
y_high = int(sys.argv[6])

cube_2 = cube_2[x_low-1:x_high,y_low-1:y_high,:,:]

cube_2 = np.transpose(cube_2,(1,0,2,3))

print cube_1.shape
print cube_2.shape

cube_1[:,:,7,:] *= np.cos(cube_1[:,:,10,:])
cube_2[:,:,7,:] *= np.cos(cube_2[:,:,10,:])
x = int(sys.argv[7])
y = int(sys.argv[8])

plt.figure(figsize=[5.5,12.0])


plt.subplot(411)
plt.plot(cube_1[x,y,0,:],cube_1[x,y,2,:],label='Fit')
plt.plot(cube_2[x,y,0,:],cube_2[x,y,2,:],label='Original')
#plt.plot(mhd_cube[i,j,0,:],mhd_cube[i,j,2,:],label='Original MHD')
plt.xlim([-5,1])
plt.ylim([3000,10000])
plt.xlabel('$\mathrm{log\,}\\tau$')
plt.ylabel('$\mathrm{T\,[K]}$')
plt.legend()

plt.subplot(412)
plt.plot(cube_1[x,y,0,:],cube_1[x,y,7,:])
plt.plot(cube_2[x,y,0,:],cube_2[x,y,7,:])
#plt.plot(mhd_cube[i,j,0,:],mhd_cube[i,j,7,:]*np.cos(mhd_cube[i,j,10,:]))
plt.xlabel('$\mathrm{log\,}\\tau$')
plt.ylabel('$\mathrm{B_{los}\,[Gauss]}$')

plt.subplot(413)
plt.plot(cube_1[x,y,0,:],cube_1[x,y,8,:]/1E5)
plt.plot(cube_2[x,y,0,:],cube_2[x,y,8,:]/1E5)
plt.xlabel('$\mathrm{log\,}\\tau$')
plt.ylabel('$v_{turb}$')

plt.subplot(414)
plt.plot(cube_1[x,y,0,:],cube_1[x,y,9,:]/1E5)
plt.plot(cube_2[x,y,0,:],cube_2[x,y,9,:]/1E5)
#plt.plot(mhd_cube[i,j,0,:],mhd_cube[i,j,9,:]/1E5)
plt.xlabel('$\mathrm{log\,}\\tau$')
plt.ylabel('$v_{los}$')

plt.tight_layout()

plt.savefig('atmos_comp',fmt='png',bbox_inches='tight')
