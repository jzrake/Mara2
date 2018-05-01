#!/usr/bin/python
def plot(primitive):
	import numpy as np
	import matplotlib.pyplot as plt
	import sys, os, h5py, re

	#using globbing and sys, but not argparse
	h5files, fileobj, primitive_arrays = [], [], []

	#get runtime parameters from h5 file
	paramfile = h5py.File('chkpt.0000.h5','r')['/user/']
	getParams = lambda param: np.array(paramfile[param])

	paramdict = {
	'gamma' : getParams('gamma'),
	'h' : getParams('h'),
	'g0' : getParams('g0'),
	'tfinal' : getParams('tfinal'),
	'N' : getParams('N')
	}

	#compute free-fall time and sound-crossing time
	ff =  10.0**(1.5)/(2.0*paramdict["g0"])**(0.5)
	sc =  10.0/(paramdict["h"]*paramdict["g0"])**0.5
	fraction = paramdict["tfinal"] / sc

	#if sys.argv is length 2 (i.e. __name__ and primitive), create txt file with all the h5 file names, get to work on those
	if len(sys.argv) == 2:
		os.system("touch myfiles.txt | ls *.h5 > myfiles.txt")
		file = open("myfiles.txt","r")
		for line in file:
			h5files.append(line[:-1])
			fileobj.append(h5py.File(line[:-1],'r'))
		os.system("rm myfiles.txt")

	else:
		h5files = sys.argv[2:]
		for string in h5files:
			fileobj.append(h5py.File(string,'r'))

	#define x_array[0:N] and primitive array to use in plt.plot()
	x_array = fileobj[0]['/mesh/points/x'][0:paramdict['N']]
	for file in fileobj:
		primitive_arrays.append(file['/primitive/'+ primitive])

	for i in range(0,len(primitive_arrays)):
		#plt.plot(x_array,primitive_arrays[i],label = "t = "+str(round(int(h5files[i][-6:-3])/441.0*fraction,2))+"sc")
		plotlabel = re.search(".*\.([0-9]*)\..*", h5files[i]).group(1)
		plt.plot(x_array,primitive_arrays[i], label = plotlabel)

	plt.title(primitive + " plot t = 0 to " + str(round(fraction,2)) + "sc" + "; ff = " + str(round(ff,2)) + "; sc =" + str(round(sc,2)))
	plt.xlabel("Radius (r)"); plt.ylabel(primitive)
	plt.savefig(primitive+".png") #save plot
	plt.legend(); plt.show()


def main():
	import sys
	primitive = sys.argv[1]
	plot(primitive)


if __name__ == "__main__":
	main()
