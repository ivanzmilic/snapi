import pyana #for writing final result
from astropy.io import fits 
import numpy as np 
import matplotlib.pyplot as plt # we will want to visualize in the meantime
import sys

output_filename = sys.argv[1]

NX = 288
NY = 288
ND = 100

spinor_cube = np.zeros([12,NX,NY,ND])

temp_input = fits.open('result_'+'0'+'.200000.fits')
pp = temp_input[0].data

h = np.linspace(-650.0,340.0,100)
h *= 1E5

# h is sorted
spinor_cube[1,:,:,:] = h


# Temperature:
temp_input = fits.open('eosT.200000.fits')
p = temp_input[0].data

for i in range(0,NX):
	for j in range(0,NY):
		spinor_cube[2,i,j,:] = p[i,:,j]

#Pressure:
temp_input = fits.open('eosP.200000.fits')
p = temp_input[0].data

for i in range(0,NX):
	for j in range(0,NY):
		spinor_cube[3,i,j,:] = p[i,:,j]

print p[0,:,0]

#density
temp_input = fits.open('result_0.200000.fits')
pp = temp_input[0].data


#velocity:
temp_input = fits.open('result_2.200000.fits')
p = temp_input[0].data
for i in range(0,NX):
	for j in range(0,NY):
		spinor_cube[9,i,j,:] = p[i,:,j]/pp[i,:,j]

#magnetic field:
temp_input = fits.open('result_5.200000.fits')
Bx = temp_input[0].data
temp_input = fits.open('result_6.200000.fits')
Bz = temp_input[0].data
temp_input = fits.open('result_7.200000.fits')
By = temp_input[0].data

for i in range(0,NX):
	for j in range(0,NY):
		spinor_cube[7,i,j,:] = np.sqrt(Bx[i,:,j]**2.0+By[i,:,j]**2.0+Bz[i,:,j]**2.0)
		spinor_cube[10,i,j,:] = np.arccos(Bz[i,:,j]/spinor_cube[7,i,j,:])
		spinor_cube[11,i,j,:] = np.arctan(By[i,:,j]/Bx[i,:,j])

d=66
plt.clf()
plt.cla()
plt.subplot(221)
plt.imshow(spinor_cube[2,:,:,d])
plt.colorbar()
plt.subplot(222)
plt.imshow(spinor_cube[7,:,:,d])
plt.colorbar()
plt.subplot(223)
plt.imshow(spinor_cube[10,:,:,d])
plt.colorbar()
plt.subplot(224)
plt.imshow(spinor_cube[9,:,:,d]/1e5)
plt.colorbar()
plt.savefig('test_0.png',fmt='png')

#finally we need to flip it and write it down:
spinor_cube = spinor_cube[:,:,:,::-1]
pyana.fzwrite(output_filename,spinor_cube,0,'placeholder')
	