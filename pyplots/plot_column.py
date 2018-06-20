import matplotlib
matplotlib.use('Agg')
import numpy as np 
import scipy.signal as sgn
import matplotlib.pyplot as plt 
import sys 
import pyana 

#print(plt.style.available)
plt.style.use('seaborn-colorblind')
#matplotlib.rcParams['font.family'] = 'serif'
#matplotlib.rcParams['font.serif'] = ['Tahoma']
matplotlib.rcParams['font.size'] = 13

nodes_in = sys.argv[1]
atmos_in = sys.argv[2]

temp = pyana.fzread(nodes_in)
nodes = temp["data"]
temp = pyana.fzread(atmos_in)
atmos = temp["data"]

i = int(sys.argv[3])
j = int(sys.argv[4])

T_nodes = [-2.5,-1.4,-0.6,0.0]
B_nodes = [-2.9,-0.5]

plt.clf()
plt.cla()
plt.subplot(211)
plt.plot(atmos[i,j,0],atmos[i,j,2])
plt.plot(T_nodes,nodes[0:4,j,i],'o',color='red',markersize=10)
plt.ylabel('$\mathrm{Temperature\,[K]}$')
plt.xlabel('$\log\\tau$')
plt.xlim([-4,1])
plt.subplot(212)
plt.plot(atmos[i,j,0],atmos[i,j,7])
plt.plot(B_nodes,nodes[7:9,j,i],'o',color='red',markersize=10)
plt.ylabel('$\mathrm{B\,[Gauss]}$')
plt.xlabel('$\log\\tau$')
plt.xlim([-4,1])
plt.ylim([0,3000])
plt.tight_layout()

plt.savefig('column_'+str(i)+'_'+str(j),fmt='png',bbox_inches='tight')
plt.savefig('column_'+str(i)+'_'+str(j)+'.eps',fmt='eps',bbox_inches='tight')