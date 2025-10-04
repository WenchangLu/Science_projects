#!/usr/bin/python3

import sys

import math
import numpy as np
import matplotlib.pyplot as plt
import os
import io

if len(sys.argv) > 1:
    with open(sys.argv[1], "r") as f:
        sts_data_lines = f.read().split("\n")
        #readlines will have \n in each list
else:
    print("need the stm cube file")  
    exit(0)

sts_data_lines = [s for s in sts_data_lines if '#' not in s]
num_points = int(sts_data_lines[0].split()[0])
pos_start = float(sts_data_lines[0].split()[1])
pos_end = float(sts_data_lines[0].split()[2])
delta_d = (pos_end - pos_start)/float(num_points-1)
E_points = int(sts_data_lines[1].split()[0])
E_min = float(sts_data_lines[1].split()[1])
E_max = float(sts_data_lines[1].split()[1])
delta_e = (E_max - E_min)/float(E_points -1)

vol_data = ''.join(sts_data_lines[2:]).split()
sts_1d = np.array(vol_data, dtype= float)

sts_2d = sts_1d.reshape(num_points,E_points)

#color_map = st.sidebar.radio("color map", ["hot", "inferno", "Greys_r"])
color_map = "hot"
color_map = "inferno"

#spline_type = col3.checkbox("spline interpolation for image", False)
X = np.ndarray([num_points, E_points], dtype = float)
Y = np.ndarray([num_points, E_points], dtype = float)
for j in range(num_points):
    for i in range(E_points):
        X[j][i] = pos_start + j * delta_d
        Y[j][i] = E_min + i * delta_e

fig, ax = plt.subplots(figsize=(15,15))

    #pos = ax.imshow(rho_ext, cmap =color_map, x_ax = X, y_ax = Y)
pos = ax.pcolormesh(X, Y, sts_2d, cmap =color_map)
#ax.set_aspect((ymax-ymin)/(xmax-xmin))
ax.set_aspect("equal")
ax.set_xlabel('X($\mathrm{\AA}$)')
ax.set_ylabel('Y($\mathrm{eV}$)')

cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])

plt.colorbar(pos, cax=cax)


fn='sts.png'
plt.savefig(fn, format='png', dpi =600)

