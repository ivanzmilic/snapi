from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import colorcet as cc

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyana
import numpy as np 
import sys
from scipy.interpolate import interp1d
from matplotlib_scalebar.scalebar import ScaleBar
from scipy.stats.stats import pearsonr
from mpl_toolkits.axes_grid1 import make_axes_locatable

node_file = sys.argv[1]
NT = int(sys.argv[2])
NVt = int(sys.argv[3])
NVs = int(sys.argv[4])
NB = int(sys.argv[5])
Ntheta = int(sys.argv[6])
Nphi = int(sys.argv[7])

output_file = sys.argv[8]

temp = pyana.fzread(node_file)
nodes = temp["data"]
dims = nodes.shape
NX = dims[1]
NY = dims[2]
NN = dims[0]

scale = 10.0

fig, axes = plt.subplots(nrows=NN,ncols=1,figsize=(scale,scale*NY/NX*NN))
image_no = 0
for i in range(0,NT):
	m = np.mean(nodes[i])
	s = np.std(nodes[i])
	ax = axes.flat[image_no]
	im = ax.imshow(nodes[i].transpose(),origin='lower',vmin = m-3*s,vmax=m+3*s,cmap='hot')
	image_no +=1

for i in range(NT,NT+NVt+NVs):
	nodes[i] /= -1E5
	m = np.mean(nodes[i])
	s = np.std(nodes[i])
	ax = axes.flat[image_no]
	im = ax.imshow(nodes[i].transpose(),origin='lower',vmin = m-3*s,vmax=m+3*s,cmap='coolwarm')
	image_no +=1

for i in range(NT+NVt+NVs, NT+NVt+NVs+NB):
	B = nodes[i] * np.cos(nodes[-1])
	m = np.mean(B)
	s = np.std(B)
	ax = axes.flat[image_no]
	im = ax.imshow(nodes[i].transpose(),origin='lower',cmap='coolwarm')
	image_no+=1


plt.tight_layout()
plt.savefig(output_file+'.eps',fmt='eps',bbox_inches='tight')
plt.savefig(output_file+'.png',fmt='png',bbox_inches='tight')