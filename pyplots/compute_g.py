import sys
s = float(sys.argv[1])
l = float(sys.argv[2])
j = float(sys.argv[3])

gls = 1.5 + (s*(s+1.0) - l*(l+1.0)) / 2.0 / j / (j+1.0)

print gls