#!/usr/bin/python3

import sys

import math
import numpy as np
import matplotlib.pyplot as plt
import os
import io

if len(sys.argv) > 2:
    with open(sys.argv[1], "r") as f:
        rho_data_1 = f.read()
        rho_data_1_lines = rho_data_1.split("\n")
    with open(sys.argv[2], "r") as f:
        rho_data_2 = f.read()
        rho_data_2_lines = rho_data_2.split("\n")
else:
    print("need two stm cube files")  
    exit(0)


rho_iso = 1.0e-6
xcells = 2
ycells = 2

rho_data_1_lines = [s for s in rho_data_1_lines if '#' not in s]
rho_data_2_lines = [s for s in rho_data_2_lines if '#' not in s]
num_atoms = int(rho_data_1_lines[0].split()[0])
real_or_complex = int(rho_data_1_lines[0].split()[4])
if real_or_complex == 2:
    print("STM data cannot be complex, check the .cube file")
    exit(0)
nx = int(rho_data_1_lines[1].split()[0])
ny = int(rho_data_1_lines[2].split()[0])
nz = int(rho_data_1_lines[3].split()[0])

atom_zpos = []
bohr = 0.529177
for i in range(num_atoms):
    atom_zpos.append(float(rho_data_1_lines[i+4].split()[4]) * bohr)


vec = np.ndarray([3,3], dtype = float)
for i in range(3):
    for j in range(3):
        vec[i][j] = float(rho_data_1_lines[1+i].split()[j+1]) * bohr

if abs(vec[0][1]) > 1.0e-5 or abs(vec[0][2]) > 1.0e-5:
    print(vec[0][1], vec[0][2])
    print("the first vector need to be along x axis")
    exit(0)

zmax = float(max(atom_zpos))
zmin = float(min(atom_zpos))
z_length = vec[2][2] * nz

z_vacuum = z_length -(zmax - zmin)
z_stm_b = float(zmax + 2.0)
z_stm_t = float(zmax + 0.5 * z_vacuum)
iz_start = int(zmax/vec[2][2])
iz_end = int( (z_stm_t)/vec[2][2])
z_stm_ave = (z_stm_b + z_stm_t) * 0.5
vol_data_1 = ''.join(rho_data_1_lines[num_atoms+4:]).split()
rho_1_1d = np.array(vol_data_1, dtype= float)

rho_1_3d = rho_1_1d.reshape(nx,ny,nz*real_or_complex)
vol_data_2 = ''.join(rho_data_2_lines[num_atoms+4:]).split()
rho_2_1d = np.array(vol_data_2, dtype= float)

rho_2_3d = rho_2_1d.reshape(nx,ny,nz*real_or_complex)
rho_xy = np.ndarray([nx,ny], dtype = float)

rho_3d = rho_1_3d - rho_2_3d

#STM_mode = st.sidebar.radio("STM mode", ["constant current", "constant height"])
STM_mode = "constant current"

if STM_mode == "constant height":
    z_height = st.sidebar.slider("choose the xy plane with z height", z_stm_b, z_stm_t, z_stm_b)
    z_plane = int(z_height/vec[2][2])

    rho_xy = rho_3d[:,:, z_plane]
else:
    rho_min = float(min(rho_3d[:,:, iz_end].reshape(nx*ny)))
    rho_max = float(max(rho_3d[:,:, iz_end].reshape(nx*ny)))
    print("rho_max min", rho_max,rho_min)
    rho_iso = rho_max*2.0

    print("isosurface: ", rho_iso)
    #print(rho_iso, iz_end)
    for ix in range(nx):
        for iy in range(ny):
           for iz in range(iz_end, 0, -1):
               if float(rho_3d[ix, iy, iz]) > rho_iso:
                  iz_tem = iz
                  break
           height = float(iz_tem) + (rho_3d[ix,iy, iz_tem] - rho_iso)/(rho_3d[ix,iy,iz_tem] - rho_3d[ix,iy,iz_tem+1])       
           rho_xy[ix][iy] = height * vec[2][2]
           #print(ix,iy,iz_tem, rho_xy[ix][iy], rho_3d[ix,iy,iz_tem], rho_3d[ix,iy,iz_tem+1])

 
#color_map = st.sidebar.radio("color map", ["hot", "inferno", "Greys_r"])
color_map = "hot"
color_map = "inferno"

#spline_type = col3.checkbox("spline interpolation for image", False)
X = np.ndarray([nx *xcells, ny *ycells], dtype = float)
Y = np.ndarray([nx *xcells, ny *ycells], dtype = float)
rho_ext = np.ndarray([nx *xcells,ny *ycells], dtype = float)
a_length = vec[0][0] * nx
for j in range(ny * ycells):
    for i in range(nx * xcells):
        X[i][j] = i * vec[0][0] + j * vec[1][0]
        Y[i][j] = i * vec[0][1] + j * vec[1][1]
        rho_ext[i][j] = float(rho_xy[i%nx][j%ny])

fig, ax = plt.subplots(figsize=(15,15))

    #pos = ax.imshow(rho_ext, cmap =color_map, x_ax = X, y_ax = Y)
pos = ax.pcolormesh(X, Y, rho_ext, cmap =color_map)
#ax.set_aspect((ymax-ymin)/(xmax-xmin))
ax.set_aspect("equal")
ax.set_xlabel('X($\mathrm{\AA}$)')
ax.set_ylabel('Y($\mathrm{\AA}$)')

cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])

plt.colorbar(pos, cax=cax)


fn=sys.argv[1].replace('cube', 'png')
plt.savefig(fn, format='png', dpi =600)

