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

# a numpy array for file numbers to read in (w/o padded zeros)
nums = np.arange(1,2,1)

datadir = './data/' # path to the directory that has the id*/ directories
probid = 'Par_Strat3d' # from athinput file
filesuff = '.ds.lis'  # suffix for the file names, from athinputfile

for num in nums:
	filenum = '{0:04d}'.format(num)

	# read in single lis file from id0/ dir
	filename_sing = datadir+'id0/'+probid+'.'+filenum+filesuff
	data, grpropertyrad = aP.read_single_lis(filename_sing)
	# this routine is particularly useful if the .lis files are already
	# joined using the join_lis.c utility available with Athena
	print('Num. particles one proc 0: ', len(data['px1']))


	# read all lis files from the id*/ directories
	data, grpropertyrad = aP.read_lis(probid, datadir, filenum, filesuff,
		numprocs)
	print('Num. particles all procs: ', len(data['px1']))


	downsample = 1 
	# if you want to downsample the number of particles for the 
	# particle plot below

	# Simple plot of the particle positions in x1-x2 plane
	plt.plot(data['px1'][::downsample], data['px2'][::downsample],
		'.', ms=1)
	plt.xlabel('particle x1 position')
	plt.ylabel('particle x2 position')
	plt.show()
