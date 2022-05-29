'''
	A script for reading in .bin data from Athena (C-version)
	using the athenaPy module.

	run with python <script name>, no arguments required.
	Changing paramters is done within the script.

	Requires numpy and matplotlib.

	J. Rucska, McMaster Univeristy, 2022.
'''
import athenaPy as aP
import numpy as np
import matplotlib.pyplot as plt

# number of mpi processes: [NGrid_x1, NGrid_x2, NGrid_x3]
numprocs = np.asarray([4,1,2])

# number of grid cells in each dimesion
Nx1 = 32 
Nx2 = 32
Nx3 = 32
Nxs = np.asarray([Nx1,Nx2,Nx3])

# a numpy array for file numbers to read in (w/o padded zeros)
nums = np.arange(1,2,1)

# Flags for in particles or self gravity is on. This changes the format of 
# the Athena binary outputs
bPar = True;
bGrav = True;

# strings for files to be read in
probid = 'Par_Strat3d' # from athinput file
filepath = './data/' # path to the directory that has the id*/ directories

for num in nums:
	filenum = '{0:04d}'.format(num)

	data = aP.read_bin(probid, filepath, filenum,
		Nxs, numprocs, bPar, bGrav)
	# output is a dictionary. Available keys are:
	#	x1, x2, x3, d (gas density),
	#	v1, v2, v3 (gas velocities).
	#	If used: dpar (dust density),
	#	m1par, m2par, m3par (dust momenta),
	#	and phi (gravitational potential, from self-grav).

	d = data['dpar']
	print(d.shape)
	plt.imshow(d[:,16,:])
	plt.show()
