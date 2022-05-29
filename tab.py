import numpy as np
import sys

# ---------------------------------------------------------------------------- #
# ----------------------------- 3D routines ---------------------------------- #
# ---------------------------------------------------------------------------- #

def read_tab_3D(probid, datapath, filenum, Nx, xlims,
	outid=None, numprocs=[1,1,1], bPrim=False, bPar=False, bGrav=False):
	"""Reads athena bin file into dictionay.
	
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
	numprocs : array_like
		array for number of MPI process used:
		[NGrid_x1, NGrid_x2, NGrid_x3].
		Default [1,1,1] for if simulation was run in serial
	bPar : bool, optional
		Boolean for if particle module was used
	bGrav : bool, optional
		Boolean for if particle module was used
	bDoublePres : bool, optional
		True if double precision was used in data output

	Returns
	-----------
	dict
		dictionary of 3D numpy arrays. Available keys are:
		x1, x2, x3, d (gas density),
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


		numproc_x1 = numprocs[0]
		numproc_x2 = numprocs[1]
		numproc_x3 = numprocs[2]
		bDoublePres = 1
	else:
		datas = np.zeros((Nx1,Nx2,Nx3))
		nVar = 1

	cnt3 = 0
	for idn3 in range(numproc_x3):	
		cnt2 = 0
		for idn2 in range(numproc_x2):
			cnt1 = 0
			for idn1 in range(numproc_x1):

				idno = idn1+numproc_x1*idn2+numproc_x1*numproc_x2*idn3

				filepath = datapath
				# need to reset filepath to string from input
				# arg. list
				afterPre = ''	
				if( numproc_x1 == 1 and
					numproc_x2 == 1 and
					numproc_x3 == 1):
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
			# Won't necessarily have a bin at every time step of the output, will have bin at 0000 always, I think.

				filename = filepath + probid + fileEnd
				filenamebin = filepath + probid + fileEndbin

				ns, ixs = parse_bin_for_inds(filenamebin,
					xlims)
				nxp = ns[0]; nyp = ns[1]; nzp = ns[2];
	
				print(filename)
				datap = parse_tab_3D(filename, ns, bPrim, nVar)

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


						#datas[cnt1:nxp+cnt1,Nx2-cnt2-jj-1,Nx3-cnt3-kk-1] = datap[(jj*nxp)+(kk*nxp*nyp):(jj+1)*nxp+kk*nxp*nyp]

				cnt1 += nxp
			cnt2 += nyp
		cnt3 += nzp


	#x1s = bin3Dflips(x1s)
	#x2s = bin3Dflips(x2s)
	#x3s = bin3Dflips(x3s)
	#dps = bin3Dflips(dps)
	#V1s = bin3Dflips(V1s)
	#V2s = bin3Dflips(V2s)
	#V3s = bin3Dflips(V3s)
	#if(bPar):
	#	dpars  = bin3Dflips(dpars)
	#	M1pars = bin3Dflips(M1pars)
	#	M2pars = bin3Dflips(M2pars)
	#	M3pars = bin3Dflips(M3pars)
	#if(bGrav):
	#	Phis = bin3Dflips(Phis)

	if(bPrim):
		# Store global 3D arrays in a single dictionary, return that
		data = {}
		data['x1'] = x1s
		#data['x2'] = x2s
		#data['x3'] = x3s
		data['d'] = dps
		data['v1'] = V1s
		data['v2'] = V2s
		data['v3'] = V3s

		if(bPar):
			data['dpar'] = dpars
			data['m1par'] = M1pars
			data['m2par'] = M2pars
			data['m3par'] = M3pars
		if(bGrav):
			data['phi'] = Phis
	else:
		data = np.copy(datas)

	return data


def parse_tab_3D(filename, ns, bPrim, nVar):
	"""test new docstring: Read in data from athena file.
	
	"""
	try:
	  file = open(filename,'r')
	except:
	  print("Couldn't open file, tried: ", filename)
	  raise SystemExit

	print(ns)

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

# ---------------------------------------------------------------------------- #
# ----------------------------- 2D routines ---------------------------------- #
# ---------------------------------------------------------------------------- #

def read_tab_2D(probid, datapath, filenum, outid, Nx, xlims, iDim,
	numprocs=[1,1,1]):
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
		print(filename.replace(datapath, ''))

		binfilename = filename.replace(outid + '.tab', 'bin')
		binfilename = binfilename.replace(filenum, '0000')

		# pull indices that belong to this processor from bin file
		ns, ixs = parse_bin_for_inds(binfilename, bDoublePres, xlims)

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
