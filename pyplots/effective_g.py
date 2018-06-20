import numpy as np 
import sys
g1 = float(sys.argv[1])
g2 = float(sys.argv[2])
j1 = float(sys.argv[3])
j2 = float(sys.argv[4])

geff = 0.5*(g2-g1) + 0.25*(g2-g1)*(j2*(j2+1.0)-j1*(j1+1.0))
print geff

