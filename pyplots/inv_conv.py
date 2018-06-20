import numpy as np 
import sys

l = float(sys.argv[1]) #this is input wavelength in A
I = float(sys.argv[2]) #this is input intensity in SI units, RH style
                       # J m^-2 s^-1 Hz^-1
c = 2.99792E10 #cm/s

l/=1E8 #convert to cm

I_out = I*c/l/l*1E3

print I_out
