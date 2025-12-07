#!/usr/bin/python3

import sys

import math
import numpy as np
import matplotlib.pyplot as plt
import os
import io

if len(sys.argv) > 1:
    with open(sys.argv[1], "r") as f:
        data_lines = f.read().split("\n")
        #readlines will have \n in each list
else:
    print("need the stm cube file")  
    exit(0)

E_range = [-4.0, 4.0]
Epoints = 801
num_grid = 448

data_lines = [s for s in data_lines if '#' not in s]
data_lines = [s for s in data_lines if s != ""]
print(len(data_lines))
X = np.array([s.split()[0] for s in data_lines], dtype = float)
Y = np.array([s.split()[1] for s in data_lines], dtype = float)
Z = np.array([s.split()[2] for s in data_lines], dtype = float)

Z = Z /2000.0
X2d = X.reshape(Epoints,num_grid)
Y2d = Y.reshape(Epoints,num_grid)
Z2d = Z.reshape(Epoints,num_grid)

#color_map = st.sidebar.radio("color map", ["hot", "inferno", "Greys_r"])
color_map = "inferno"
color_map = "hot"
color_map = "bwr"

fig, ax = plt.subplots(figsize=(15,15))

#pos = ax.pcolormesh(X2d, Y2d, Z2d, cmap =color_map)
pos = ax.pcolormesh(X2d, Y2d, Z2d, cmap =color_map, vmax = 30)

#ax.set_aspect((pos_end-pos_start)/(E_max-E_min))
#ax.set_aspect((E_range[1]-E_range[0])/(pos_end-pos_start))
#ax.set_aspect("equal")
#plt.xticks(fontsize=12) 
#plt.yticks(fontsize=12) 

# Or using the Axes object
ax.tick_params(axis='both', which='major', labelsize=24)
ax.set_ylabel('E(eV)', fontsize = 24)
ax.set_xlabel('Z($\AA$)', fontsize = 24)
#plt.xlim(20,40)
plt.ylim(-2,2)

cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])

cbar = fig.colorbar(pos, cax=cax)
cbar.ax.tick_params(labelsize=24) 
cbar.set_label("DOS(arbitrary units)", fontsize = 24)


fn='pdos.png'
plt.savefig(fn, format='png', dpi =600)

