import numpy as np

ll1 = np.linspace(7698.0, 7700.0, 201)
ll2 = np.linspace(8540.0, 8544.0, 401)
ll3 = np.linspace(10826.0, 10828.0, 201)

ll = np.concatenate([ll1, ll2, ll3])

np.savetxt("lgrid_k_ca_si.dat", ll, fmt='%12.6e', header=str(len(ll)), comments='')