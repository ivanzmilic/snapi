import numpy as np 
import sys

lambda_air = float(sys.argv[1])

h = 6.62607004E-27
c = 2.99792458E10
ev = 1.60218e-12

s = 10000.0/lambda_air
n = 1.0 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s**2.0) + 0.0001599740894897 / (38.92568793293 - s**2.0)
lambda_vac = lambda_air * n
print lambda_vac
print h*c/lambda_vac / ev * 1E8