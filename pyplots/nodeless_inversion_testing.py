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

U,w,V = np.linalg.svd(H)

J_T = J[:,0:ND]
J_V = J[:,ND:2*ND]

H_T = np.dot(J_T.transpose(),J_T)
H_V = np.dot(J_V.transpose(),J_V)
for i in range(0,ND):
	H_T[i,i] *= 1.001
	H_V[i,i] *= 1.001

#small = np.where(w<1E-3*max(w))
#print small

UT,wT,VT = np.linalg.svd(H_T)
UV,wV,VV = np.linalg.svd(H_V)

plt.clf()
plt.cla()

for i in range(0,5):
	#w_inv = 1.0/w
	#w_inv[i:] = 0.0
	#w_inv[small] = 0.0;
	#W = np.diag(w_inv)

	#H_inv = np.dot(V.transpose(),np.dot(W,U.transpose()))
	#b = np.dot(JT,res[2])
	#correction = np.dot(H_inv,b)

	plt.subplot(121)
	w_inv = 1.0/wT
	w_inv[i:] = 0.0
	W = np.diag(w_inv)
	H_inv = np.dot(VT.transpose(),np.dot(W,UT.transpose()))
	b = np.dot(J_T.transpose(),res[2])
	correction = np.dot(H_inv,b)
	plt.plot(correction,label=str(i))
	plt.ylim([-2000.0,2000.0])
	
	print 'TEMP'
	for j in range(0,ND):
		print j,correction[j]
	
	plt.subplot(122)
	w_inv = 1.0/wV
	w_inv[i:] = 0.0
	W = np.diag(w_inv)
	H_inv = np.dot(VV.transpose(),np.dot(W,UV.transpose()))
	b = np.dot(J_V.transpose(),res[2])
	correction = np.dot(H_inv,b)	
	plt.plot(correction/1E5,label=str(i))
	plt.legend()

	print 'VEL'
	for j in range(0,ND):
		print j,correction[j]/1E5
	

plt.savefig('corrections.png')

#print correction

#print H.shape
#print b.shape


plt.clf()
plt.cla()
plt.imshow(rf[0,:,:,0])
plt.savefig('justesting',fmt='png')