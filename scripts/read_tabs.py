'''
	A script for reading in .tab data from Athena (C-version)
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

# domain limits: [x1min, x2min, x3min]
xmins = [-0.1, -0.1, -0.1]

# a numpy array for file numbers to read in (w/o padded zeros)
nums = np.arange(1,2,1)

# strings for files to be read in
probid = 'Par_Strat3d' # from athinput file
datapath = './data/' # path to the directory that has the id*/ directories
outid = 'dparxy' # from athinput file

iDim = 3 # what dimension 2D data is collapsed in
# using x3 = : in <output> block in inputfile means to use iDim = 3

for num in nums:
	filenum = '{0:04d}'.format(num)

	data = aP.read_tab_3D(probid, datapath, filenum, Nxs, xmins,
		outid=None, numprocs=numprocs, 
		bPrim=True, bPar=True, bGrav=True)
	# Reading in prim data (bPrim = True), output is a dictionary.
	#	Available keys are:
	#	x1, d (gas density), (x2, x3 not implemented currently)
	#	v1, v2, v3 (gas velocities).
	#	If used: dpar (dust density),
	#	m1par, m2par, m3par (dust momenta),
	#	and phi (gravitational potential, from self-grav).

	data_2D = aP.read_tab_2D(probid, datapath, filenum, outid, Nxs,
		xmins, iDim, numprocs)

	print('data shapes: ', data['dpar'].shape, data_2D.shape)

	fig, axs = plt.subplots(2,1)
	axs[0].imshow(data['d'][:,16,:])
	axs[1].imshow(data_2D)
	plt.show()

