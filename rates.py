import numpy as np 

np.set_printoptions(precision = 15)

rates = np.loadtxt("rates.txt", unpack = True)

rates = rates.transpose()
rates_0 = rates[0:3,0:3]
rates_p = rates[3:6,0:3]
rates_m = rates[6:9,0:3]

b_0 = rates[0:3,3]
b_p = rates[3:6,3]
b_m = rates[6:9,3]

pops = np.loadtxt("pops.txt", unpack = True)
pops = pops.transpose()
p_0 = pops[0,:]
p_p = pops[1,:]
p_m = pops[2,:]

#print p_0
#print np.linalg.solve(rates_0, b_0)

#print p_p
#print np.linalg.solve(rates_p, b_p)

#print p_m
#print np.linalg.solve(rates_m, b_m)

d_rates = rates_p - rates_m

d_b = b_p - b_m
#print d_b


beta = np.zeros(3)
for i in range(0,3):
	for j in range(0,3):
		beta[i] += -d_rates[i,j] * p_0[j]

beta += d_b

print np.linalg.solve(rates_0, beta)
print p_p - p_m





