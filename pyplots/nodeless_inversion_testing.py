import numpy as np 
import matplotlib.pyplot as plt 
import sys 

rf_file = sys.argv[1]
res_file = sys.argv[2]

rf = np.loadtxt(rf_file)

ND = 57
rf = rf[:,2:6].reshape(7,57,194,4)

res = np.loadtxt(res_file,unpack=True)
dims2 = res.shape
NL = dims2[1]

J = np.zeros([NL,2*ND])

for d in range(0,ND):
	J[:,d] = rf[0,d,:,0]
	J[:,d+ND] = rf[3,d,:,0]

JT = J.transpose()
H = np.dot(JT,J)
for i in range(0,2*ND):
	H[i,i] *= 1.001
	
b = np.dot(JT,res[2])

correction = np.linalg.solve(H,b)

print correction

print H.shape
print b.shape


plt.clf()
plt.cla()
plt.imshow(rf[0,:,:,0])
plt.savefig('justesting',fmt='png')