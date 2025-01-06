#!/usr/bin/python3

import sys

import math
import numpy as np
import matplotlib.pyplot as plt
import os
import io

dir_list = ["bin30"]

with open(dir_list[0] + "/INCAR", "r") as f:
        incar_lines = f.readlines()

for line in incar_lines:
    if "EINT" in line:
        E_range = line.split("=")[1].split()

try:
    Energy = (float(E_range[0]) + float(E_range[1]))*0.5
except:
    print("no EINT in INCAR")
    exit(0)

with open(dir_list[0] + "/PARCHG", "r") as f:
    parchg_lines = f.readlines()

scale = float(parchg_lines[1])
latt = np.array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
latt[0] = [scale *float(a) for a in parchg_lines[2].split()]
latt[1] = [scale *float(a) for a in parchg_lines[3].split()]
latt[2] = [scale *float(a) for a in parchg_lines[4].split()]
num_atoms = 0
for a in parchg_lines[6].split():
    num_atoms += int(a)

grids = [int(a) for a in parchg_lines[num_atoms + 9].split()]
nx = grids[0]
ny = grids[1]
nz = grids[2]
rho_list = []
del parchg_lines[0:num_atoms+10]
fname = dir_list[0] + "/PARCHG"
rho = np.loadtxt(fname, dtype = float, skiprows=num_atoms+10)
if rho.size != nx*ny*nz:
    print(rho.sizei, nx, ny, nz)
    exit(0)
rho_3d = rho.reshape((nz,ny,nx))    
i_ave = 18
delt_dx = (latt[0][0] + latt[1][0])/float(nx)
delt_dy = (latt[0][1] + latt[1][1])/float(ny)
delt_d = math.sqrt(delt_dx* delt_dx + delt_dy*delt_dy)
for i in range(i_ave, nx-i_ave):
    tem = 0.0
    for j in range(-i_ave, i_ave):
        for k in range(-i_ave, i_ave):
            tem += rho_3d[72][i+j][i+k]
    print(delt_d * float(i), tem)        

color_map = "inferno"

#spline_type = col3.checkbox("spline interpolation for image", False)
X = np.ndarray([nx, ny], dtype = float)
Y = np.ndarray([nx, ny], dtype = float)

vec = latt/float(nx)
for j in range(ny):
    for i in range(nx):
        X[i][j] = i * vec[0][0] + j * vec[1][0]
        Y[i][j] = i * vec[0][1] + j * vec[1][1]

fig, ax = plt.subplots(figsize=(15,15))

pos = ax.pcolormesh(X, Y, rho_3d[72], cmap =color_map)
#ax.set_aspect((ymax-ymin)/(xmax-xmin))
ax.set_aspect("equal")
ax.set_xlabel('X($\mathrm{\AA}$)')
ax.set_ylabel('Y($\mathrm{\AA}$)')

cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])

plt.colorbar(pos, cax=cax)


plt.savefig("a.png", format='png', dpi =600)




