#!/usr/bin/python3

import sys
import os
import math
import numpy
import matplotlib.pyplot as plt

if len(sys.argv) >1:
    f = open(sys.argv[1], 'r')
else:
    f = open('dipole.dat_input.00','r')
    
damping_exp = 1.0/100.0
#damping_exp = 0.0/100.0
n_pad = 0
all_lines = f.readlines()
f.close

line = all_lines[1].split()
efield_x = abs(float(line[2]))
efield_y = abs(float(line[3]))
efield_z = abs(float(line[4]))
  
if efield_x > 0.0: 
    if efield_y > 0.0: print ('warning: both x and y field >0.0')
    if efield_z > 0.0: print ('warning: both x and z field >0.0')
    col_pick = 1
elif efield_y > 0.0: 
    if efield_z > 0.0: print ('warning: both y and z field >0.0')
    col_pick = 2
elif efield_z > 0.0:
    col_pick = 3
else:
    print ('all x, y, z field == 0')
  
      
n_data = len(all_lines) -3
fout = open('spectra.dat', 'w')
for col_pick in range(1,4):
    dipole = numpy.zeros(n_data)
    t = numpy.zeros(n_data)
    for i in range(n_data):
        line = all_lines[i + 3].split()
        t[i] = float(line[0])
        dipole[i] = float(line[col_pick]) 
        print( i, dipole[i])
    print("&")
      
    mean_dipole =  numpy.mean(dipole)

    for i in range(n_data):
        dipole[i] = (dipole[i] - mean_dipole) 
        dipole[i] = dipole[i] * math.exp(-abs(t[i]-t[0]) * damping_exp)
        print( i, dipole[i])
        

    fw = numpy.fft.rfft(dipole, n_data+n_pad)
      
    num_freq = (n_data+n_pad)//2 +1
    dt = t[1] - t[0]
    period = (n_data+n_pad) * dt
    dw = 2.0 * 3.1415926/period * 27.2114

    wmin = 0.0
    wmax = num_freq * dw
    x = numpy.linspace(wmin, wmax, num_freq)



    #  rotate spectrum as in NWCHEM
    for i in range(len(fw)):
        re = fw[i].real
        im = fw[i].imag
        r = math.sqrt (re**2 + im**2)
        angle = abs (math.atan2 (im, re))
        if angle > 3.1415927:
          print(angle)
          raise Exception ("atan2 out of range")
        fw[i] = r * math.cos(angle) + r*math.sin(angle) * 1j
        if x[i] < 5.0:
            fout.write("%f  %f\n"%(x[i], fw[i].imag))
    fout.write("&\n")
        
