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

E_range = [-1.0, 1.0]
num_points_ave = 16

sts_data_lines = [s for s in sts_data_lines if '#' not in s]
num_points = int(sts_data_lines[0].split()[0])
pos_start = float(sts_data_lines[0].split()[1])
pos_end = float(sts_data_lines[0].split()[2])
delta_d = (pos_end - pos_start)/float(num_points-1)
print(delta_d)
E_points = int(sts_data_lines[1].split()[0])
E_min = float(sts_data_lines[1].split()[1])
E_max = float(sts_data_lines[1].split()[2])
delta_e = (E_max - E_min)/float(E_points -1)

vol_data = ''.join(sts_data_lines[2:]).split()
sts_1d = np.array(vol_data, dtype= float)

sts_2d = sts_1d.reshape(num_points, E_points)

sts_ave = np.ndarray([num_points - num_points_ave, E_points], dtype = float)
for i in range(num_points - num_points_ave):
    sts_ave[i] = sum(sts_2d[i:i+ num_points_ave])

#color_map = st.sidebar.radio("color map", ["hot", "inferno", "Greys_r"])
color_map = "inferno"
color_map = "hot"
color_map = "bwr"

#spline_type = col3.checkbox("spline interpolation for image", False)
X = np.ndarray([num_points-num_points_ave, E_points], dtype = float)
Y = np.ndarray([num_points-num_points_ave, E_points], dtype = float)
for j in range(num_points-num_points_ave):
#    print("&&")
    for i in range(E_points):
        ene = E_min + i * delta_e
#        if(ene >-2.0 and ene < 2.0 and j%16 == 0):
#            print(E_min + i *delta_e, sts_2d[j][i])
        X[j][i] = pos_start + (j+num_points_ave//2) * delta_d 
        Y[j][i] = E_min + i * delta_e

fig, ax = plt.subplots(figsize=(15,15))

    #pos = ax.imshow(rho_ext, cmap =color_map, x_ax = X, y_ax = Y)
#pos = ax.pcolormesh(Y.T, X.T, sts_ave.T, cmap =color_map, vmax = 2.0)
pos = ax.pcolormesh(Y.T, X.T, np.log(sts_ave.T), cmap =color_map)
#ax.set_aspect((pos_end-pos_start)/(E_max-E_min))
ax.set_aspect((E_range[1]-E_range[0])/(pos_end-pos_start))
#ax.set_aspect("equal")
ax.set_xlabel('X(eV)')
ax.set_ylabel('Y($\AA$)')
plt.xlim(E_range[0], E_range[1])

cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])

plt.colorbar(pos, cax=cax)


fn='sts.png'
plt.savefig(fn, format='png', dpi =600)

