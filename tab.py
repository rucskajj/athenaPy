import numpy as np
import sys

# ---------------------------------------------------------------------------- #
# ----------------------------- 3D routines ---------------------------------- #
# ---------------------------------------------------------------------------- #

def read_tab_3D(probid, datapath, filenum, Nx, xlims,
	outid, numprocs=[1,1,1], 
	bPrim=False, bPar=False, bGrav=False, bDoublePres=True):
	"""Reads athena tab file into dictionary (for prim) or 3D array.
	
	Note, this routine requires at least one .bin file output from
	the simulation.

	Parameters
	-------------
	probid : str
		problem_id from athena <job> block in input file
	datapath : str
		path to data files (or id*/ directories if using MPI)
	filenum : str
		zero-padded four-digit file number
	Nx : array_like
		Number of grid points: [Nx1, Nx2, Nx3]
	xlims : array_like
		lower limits of grid domain: [x1min, x2min, x3min]
	outid: str, optional
		"id" from <output*> block in athena. If bPrim=True,
		this must be "None".
	numprocs : array_like
		array for number of MPI process used:
		[NGrid_x1, NGrid_x2, NGrid_x3].
		Default [1,1,1] for if simulation was run in serial
	bPrim : bool, optional
		True if the .tab file is for all "prim" variables
	bPar : bool, optional
		True if particle module was used
	bGrav : bool, optional
		True if self-gravity module was used
	bDoublePres : bool, optional
		True if double precision was used in bin data output

	Returns
	-----------
	dict, array
		if bPrim=False, return 3D array of desired data
		from specified .tab file.
		if bPrim=True, return dictionary of 3D numpy arrays.
		Available keys are:
		x1, d (gas density) (Note: x2, x3 not implemented
		currently),
		v1, v2, v3 (gas velocities).
		If used: dpar (dust density),
		m1par, m2par, m3par (dust momenta),
		and phi (gravitational potential, from self-grav).
	"""


	Nx1 = Nx[0]; Nx2 = Nx[1]; Nx3 = Nx[2];

	i = np.zeros((Nx1,Nx2,Nx3))
	j = np.zeros((Nx1,Nx2,Nx3))
	k = np.zeros((Nx1,Nx2,Nx3))

	# 3D array of zeroes to store data in
	if(bPrim):
		x1s = np.zeros((Nx1,Nx2,Nx3))
		x2s = np.zeros((Nx1,Nx2,Nx3))
		x3s = np.zeros((Nx1,Nx2,Nx3))

		dps = np.zeros((Nx1,Nx2,Nx3))
		V1s = np.zeros((Nx1,Nx2,Nx3))
		V2s = np.zeros((Nx1,Nx2,Nx3))
		V3s = np.zeros((Nx1,Nx2,Nx3))
		
		nVar = 7
		if(bPar):
			dpars = np.zeros((Nx1,Nx2,Nx3))
			M1pars = np.zeros((Nx1,Nx2,Nx3))
			M2pars = np.zeros((Nx1,Nx2,Nx3))
			M3pars = np.zeros((Nx1,Nx2,Nx3))
			nVar += 4
		if(bGrav):
			Phis = np.zeros((Nx1,Nx2,Nx3))
			nVar += 1
	else:
		datas = np.zeros((Nx1,Nx2,Nx3))
		nVar = 1

	# nested loop to go through all id*/ directories
	cnt3 = 0
	for idn3 in range(numprocs[2]):	
		cnt2 = 0
		for idn2 in range(numprocs[1]):
			cnt1 = 0
			for idn1 in range(numprocs[0]):

				# id*/ dir number
				idno = idn1 +\
					numprocs[0]*idn2+\
					numprocs[1]*numprocs[0]*idn3

				# construct full file name
				filepath = datapath
				afterPre = ''	
				if( numprocs[0] == 1 and
					numprocs[1] == 1 and
					numprocs[2] == 1):
					filepath += ''
					afterPre += ''
				else:
					if(idno == 0):
						filepath += 'id' + str(idno)
						afterPre += ''
					else:
						filepath += 'id' + str(idno)
						afterPre += '-id' + str(idno)
				filepath += '/'
				afterPre += '.'

				if(outid is not None):
					fileEnd = afterPre + filenum +\
						'.' + outid + '.tab'
				else:
					fileEnd = afterPre + filenum +'.tab'

				fileEndbin = afterPre + '0000' + '.bin' 

				filename = filepath + probid + fileEnd
				filenamebin = filepath + probid + fileEndbin

				# this routine requies a bin file to get
				# number of grid points on this process
				ns, ixs = parse_bin_for_inds(filenamebin,
					xlims, bDoublePres)
				nxp = ns[0]; nyp = ns[1]; nzp = ns[2];
	
				print(filename)
				# read target .tab file in this id*/ dir
				datap = parse_tab_3D(filename, ns, bPrim, nVar)

				# loop through datap from single .tab file
				# store data in global 3D array(s)
				for kk in range(nzp):		
					for jj in range(nyp):
						# integers for 3D array indices
						d1l = cnt1; d1u = cnt1+nxp;
						d2i = Nx2-cnt2-jj-1
						d3i = Nx3-cnt3-kk-1
						vl = (jj*nxp)+(kk*nxp*nyp)
						vu = (jj+1)*nxp+kk*nxp*nyp


						if(bPrim==False):
							datas[d1l:d1u,d2i,d3i]=\
							datap[vl:vu]
						else:
							x1s[d1l:d1u,d2i,d3i]=\
							datap['x1'][vl:vu]

							dps[d1l:d1u,d2i,d3i]=\
							datap['d'][vl:vu]

							V1s[d1l:d1u,d2i,d3i]=\
							datap['v1'][vl:vu]

							V2s[d1l:d1u,d2i,d3i]=\
							datap['v2'][vl:vu]

							V3s[d1l:d1u,d2i,d3i]=\
							datap['v3'][vl:vu]

							if(bPar):
								dpars[d1l:d1u,
								d2i,d3i]\
								=datap['dpar']\
									[vl:vu]

								M1pars[d1l:d1u,
								d2i,d3i]\
								=datap['m1par']\
									[vl:vu]
								M2pars[d1l:d1u,
								d2i,d3i]\
								=datap['m2par']\
									[vl:vu]
								M3pars[d1l:d1u,\
								d2i,d3i]\
								=datap['m3par']\
									[vl:vu]
							if(bGrav):
								Phis[d1l:d1u,
								d2i,d3i]\
								=datap['phi']\
									[vl:vu]
				#if(bPrim):
				#	for ii in range(nxp):
				#		for kk in range(nzp): 
				#			d1i = ii+cnt1
				#			d2l = Nx2-cnt2-nyp
				#			d2u = Nx2-cnt2
				#			d3i = Nx3-cnt3-kk-1

				#			x2s[d1i,d2l:d2u,d3i] =\
				#			np.fliplr(\
				#			[datap['x2']])[0]			
				#	for ii in range(nxp):
				#		for jj in range(nyp):
				#			d1i = ii+cnt1
				#			d2i = Nx2-cnt2-jj-1
				#			d3l = Nx3-cnt3-nzp
				#			d3u = Nx3-cnt3

				#			x3s[d1i,d2i,d3l:d3u] =\
				#			np.fliplr(\
				#			[datap['x3']])[0]

				#advance positions in global 3D array indices
				cnt1 += nxp
			cnt2 += nyp
		cnt3 += nzp


	if(bPrim):
		# Store global 3D arrays in a single dictionary, return that
		# flip 3D arrays for orientation in matplotlib imshow
		data = {}
		data['x1'] = flip_3D_array(x1s)
		#data['x2'] = flip_3D_array(x2s) currently in dev.
		#data['x3'] = flip_3D_array(x3s)
		data['d' ] = flip_3D_array(dps)
		data['v1'] = flip_3D_array(V1s)
		data['v2'] = flip_3D_array(V2s)
		data['v3'] = flip_3D_array(V3s)

		if(bPar):
			data['dpar'] = flip_3D_array(dpars)
			data['m1par'] = flip_3D_array(M1pars)
			data['m2par'] = flip_3D_array(M2pars)
			data['m3par'] = flip_3D_array(M3pars)
		if(bGrav):
			data['phi'] = flip_3D_array(Phis)
	else:
		data = np.copy(datas)

	return data

# ---------------------------------------------------------------------------- #
# ----------------------------- 2D routines ---------------------------------- #
# ---------------------------------------------------------------------------- #

def read_tab_2D(probid, datapath, filenum, outid, Nx, xlims, iDim,
	numprocs=[1,1,1], bDoublePres=True):
	"""Reads 2D slice projection outputs from athena simulations.
	
	Note, this routine requires at least one .bin file output from
	the simulation.

	Parameters
	-------------
	probid : str
		problem_id from athena <job> block in input file
	datapath : str
		path to data files (or id*/ directories if using MPI)
	filenum : str
		zero-padded four-digit file number
	outid : str
		output id from athena <output> block input block
	Nx : array_like
		Number of grid points: [Nx1, Nx2, Nx3]
	xlims : array_like
		lower limits of grid domain: [x1min, x2min, x3min]
	iDim: int
		what dimesion the slice or projection is along, values
		of 1, 2, 3 to match athena's x1, x2, x3 designation
	numprocs : array_like, optional
		array for number of MPI process used:
		[NGrid_x1, NGrid_x2, NGrid_x3]
		use default [1,1,1] if simulation was run in serial
	bDoublePres : bool, optional
		True if double precision was used in bin data output

	Returns
	-----------
	array
		2D slice of simulation data as numpy array
	"""


	Nx1 = Nx[0]; Nx2 = Nx[1]; Nx3 = Nx[2];
	
	numproc_x1 = numprocs[0]
	numproc_x2 = numprocs[1]
	numproc_x3 = numprocs[2]
	bDoublePres = True

	filelist = []

	bFound = False
	for idn3 in range(numprocs[2]):	
		for idn2 in range(numprocs[1]):
			for idn1 in range(numprocs[0]):

				# integer to count id*/ directories
				idno = idn1 + \
					numprocs[0]*idn2 +\
					numprocs[0]*numprocs[1]*idn3

				# need to reset filepath varaibles
				filepath = datapath
				probiddir = ''
		
				# construct string for path to data file
				if( numproc_x1 == 1 and
					numproc_x2 == 1 and
					numproc_x3 == 1):
					filepath += ''
					probiddir += '.'
				else:
					if(idno == 0):
						filepath += 'id0/'
						probiddir += '.'
					else:
						filepath += 'id' + \
							str(idno) + '/'
						probiddir += '-id' + \
							str(idno) + '.'

				fileEnd = probiddir + filenum + '.' + \
					outid + '.tab'
				filename = filepath + probid + fileEnd

				try:
					file = open(filename,'rb')
				except:
					continue

				# if file exists in this loops' id*/ directory,
				# add it to list of files to be opened
				filelist.append(filename)
				bFound = True
				file.close()
				# N.B. depending on the projection and grid
				# configuration, files may exist in some id*/
				# directories but not others

	if(bFound == False):
		print("Couldn't open file in any id*/ dir, last tried: \n",
			filename)
		raise SystemExit	

	inds = [1,2,3]
	inds.remove(iDim) # remove dimension data is integrated along
	inds[:] = [i-1 for i in inds] # shift 1, 2, 3's to 0, 1, 2
	outdata = np.zeros([Nx[inds[0]], Nx[inds[1]]])

	# loop through files found in main loop
	for filename in filelist:
		print(filename)

		binfilename = filename.replace(outid + '.tab', 'bin')
		binfilename = binfilename.replace(filenum, '0000')

		# pull indices that belong to this processor from bin file
		ns, ixs = parse_bin_for_inds(binfilename, xlims, bDoublePres)

		# find indices for these data in global 2D slice
		# -1 in min index calculations for indexing numpy array
		minx = ixs[inds[0]].min()-1; maxx = ixs[inds[0]].max();
		miny = ixs[inds[1]].min()-1; maxy = ixs[inds[1]].max();

		outdatabuf = np.zeros([maxx-minx, maxy-miny])
		outdatabuf = parse_single_tab_2D(filename, outdatabuf)

		# store buffer from this file in the global array
		outdata[minx:maxx, miny:maxy] = outdatabuf

	return outdata

# -------------------------------------------------------------------------- #
# ------------------------- Internal functions ----------------------------- #
# -------------------------------------------------------------------------- #

def flip_3D_array(arr):
	"""Flip 3D arrays to orient with matplotlib axes.
	
	"""
	return np.flip(np.flip(arr,axis=1),axis=2)

def parse_tab_3D(filename, ns, bPrim, nVar):
	"""Read data from single .tab file intro array or dictionary.
	
	"""
	try:
	  file = open(filename,'r')
	except:
	  print("Couldn't open file, tried: ", filename)
	  raise SystemExit

	if(bPrim): 
		headercount = 0
		datal = np.zeros([ns[0]*ns[1]*ns[2] , nVar])
	else:
		datal = np.zeros(ns[0]*ns[1]*ns[2])

	ii = 0
	for line in file:
		# prim output has 8-line header in .tab
		if(bPrim and headercount < 8):
			headercount += 1
			continue
	
		dataline = np.asarray((line.strip()).split()).astype('float64')
		datal[ii, 0:nVar] = dataline[3:3+nVar]
		ii += 1 

	if(bPrim):
		data = {}
		data['x1'] = datal[:,0]
		data['x2'] = datal[:,1]
		data['x3'] = datal[:,2]
		data['d']  = datal[:,3]
		data['v1'] = datal[:,4]
		data['v2'] = datal[:,5]
		data['v3'] = datal[:,6]

		if(nVar == 11 or nVar == 12): # bPar=True; particles are on
			data['dpar']  = datal[:,8]
			data['m1par'] = datal[:,9]
			data['m2par'] = datal[:,10]
			data['m3par'] = datal[:,11]
		if(nVar == 8 or nVar == 12): # bGrav=True; self-gravity is on
			data['phi'] = datal[:,7]
	else:
		data = datal

	return data


def parse_single_tab_2D(filename, d):
	"""Reads single .tab file into 2D numpy array.
	
	"""

	file = open(filename,'rb')

	for line in file:
		# read in data line by line
		dataline = np.asarray((line.strip()).split()).astype('float64')

		i = dataline[0].astype('int')
		j = dataline[1].astype('int')
		d[i,j] = dataline[2]
	file.close()

	return d

def parse_bin_for_inds(filename, xlims, bDoublePres=True):
	"""Reads header of single .bin file header for grid indices.
	
	"""


	try:
		file = open(filename,'rb')
	except:
		print("Couldn't find 0000.bin file, tried: \n", filename)
		raise SystemExit

	floatSize = np.float32
	filesize_float = 4
	if(bDoublePres):
		floatSize = np.float64
		filesize_float = 8

	file.seek(0,0) # start at top of file
	file.seek(4, 1) # skip coordsys

	nx,ny,nz = np.fromfile(file,dtype=np.int32,count=3)[0:3]

	nx = nx.astype(int)
	ny = ny.astype(int)
	nz = nz.astype(int)

	file.seek(4*4,1) # skip NVAR, NSCALARS, iSelfGravity, iParticles
	file.seek(filesize_float*4,1) # skip gamma1, cs, t, dt

	x1 = np.fromfile(file,dtype=floatSize,count=nx)
	x2 = np.fromfile(file,dtype=floatSize,count=ny)
	x3 = np.fromfile(file,dtype=floatSize,count=nz)

	dx1 = x1[1]-x1[0]
	# xs contains the positions of domain boundary. e.g. x[0][0] is x1min
	ix1 = np.round((x1+0.5*dx1 - xlims[0])/dx1).astype('int')
	# ix1 is now array of global/domain wide cell index

	dx2 = x2[1]-x2[0]
	ix2 = np.round((x2+0.5*dx2 - xlims[1])/dx2).astype('int')
	dx3 = x3[1]-x3[0]
	ix3 = np.round((x3+0.5*dx3 - xlims[2])/dx3).astype('int')

	file.close()

	return [nx, ny, nz], [ix1, ix2, ix3]
