#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import h5py
import argparse, sys, os

#using globbing and sys, but not argparse
h5files, fileobj, density_arrays = [], [], []
plot_count = 0

#get runtime parameters from h5 file
gamma = np.array(h5py.File('chkpt.0000.h5','r')['/user/gamma'])
h = np.array(h5py.File('chkpt.0000.h5','r')['/user/h'])
g0 = np.array(h5py.File('chkpt.0000.h5','r')['/user/g0'])
tfinal = np.array(h5py.File('chkpt.0000.h5','r')['/user/tfinal'])
N = np.array(h5py.File('chkpt.0000.h5','r')['/user/N'])

#compute free-fall time and sound-crossing time
ff =  10.0**(1.5)/(2.0*g0)**(0.5)
sc =  10.0/(h*g0)**0.5
fraction = tfinal / sc

#if sys.argv is length 1, create txt file with all the h5 file names, get to work on those
if len(sys.argv) == 1:
	os.system("touch myfiles.txt | ls *.h5 > myfiles.txt")
	file = open("myfiles.txt","r")
	for line in file:
		h5files.append(line[:-1])
		fileobj.append(h5py.File(line[:-1],'r'))

else:
	h5files = sys.argv[1:]
	for string in h5files:
		fileobj.append(h5py.File(string,'r'))

#define x_array[0:100] and density array to use in plt.plot()
x_array = fileobj[0]['/mesh/points/x'][0:N]
for file in fileobj:
	density_arrays.append(file['/primitive/density'])

for i in range(0,len(density_arrays)):
	#plt.plot(x_array,density_arrays[i],label = "t = "+str(round(int(h5files[i][-6:-3])/441.0*fraction,2))+"sc")
	plt.plot(x_array,density_arrays[i], label = str(h5files[i][6:-3]))
	plot_count += 1

plt.legend()
plt.title("Density plot t = 0 to " + str(round(fraction,2)) + "sc" + "; plot count = " + str(plot_count) + "; ff = " + str(round(ff,2)) + "; sc =" + str(round(sc,2)))
plt.xlabel("Radius (r)")
plt.ylabel("Density")
plt.show()

os.system("rm myfiles.txt")
