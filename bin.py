import numpy as np

def read_bin(probid, datapath, filenum, Nx, numprocs=[1,1,1],
	bPar=False, bGrav=False, bDoublePres=True):
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

	# 3D array of zeroes to store data in
	x1s = np.zeros((Nx1,Nx2,Nx3))
	x2s = np.zeros((Nx1,Nx2,Nx3))
	x3s = np.zeros((Nx1,Nx2,Nx3))

	dps = np.zeros((Nx1,Nx2,Nx3))
	V1s = np.zeros((Nx1,Nx2,Nx3))
	V2s = np.zeros((Nx1,Nx2,Nx3))
	V3s = np.zeros((Nx1,Nx2,Nx3))
	
	if(bPar):
		dpars = np.zeros((Nx1,Nx2,Nx3))
		M1pars = np.zeros((Nx1,Nx2,Nx3))
		M2pars = np.zeros((Nx1,Nx2,Nx3))
		M3pars = np.zeros((Nx1,Nx2,Nx3))
	
	if(bGrav):
		Phis = np.zeros((Nx1,Nx2,Nx3))

	# nested loop to go through all id*/ directories
	cnt3 = 0
	for idn3 in range(numprocs[2]):	
		cnt2 = 0
		for idn2 in range(numprocs[1]):
			cnt1 = 0
			for idn1 in range(numprocs[0]):

				# id*/ dir number
				idno = idn1 +\
					numprocs[0]* idn2 +\
					numprocs[0]*numprocs[1]*idn3

				# construct full filename from inputs
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

				fileEnd = afterPre + filenum + '.bin'
				filename = filepath + probid + fileEnd

				# read in .bin file in this id*/ dir
				nxp, nyp, nzp, datap = parse_bin(\
					filename, bDoublePres)	

				# loop through datap from single .bin file,
				# store data in global 3D arrays
				for kk in range(nzp):		
					for jj in range(nyp):
						# integers for 3D array indices
						d1l = cnt1; d1u = cnt1+nxp;
						d2i = Nx2-cnt2-jj-1
						d3i = Nx3-cnt3-kk-1
						vl = (jj*nxp)+(kk*nxp*nyp)
						vu = (jj+1)*nxp+kk*nxp*nyp

						x1s[d1l:d1u,d2i,d3i]=\
							datap['x1']
						dps[d1l:d1u,d2i,d3i]=\
							datap['d'][vl:vu]
						V1s[d1l:d1u,d2i,d3i]=\
							datap['v1'][vl:vu]
						V2s[d1l:d1u,d2i,d3i]=\
							datap['v2'][vl:vu]
						V3s[d1l:d1u,d2i,d3i]=\
							datap['v3'][vl:vu]

						if(bPar):
							dpars[d1l:d1u,d2i,d3i]\
							= datap['dpar'][vl:vu]

							M1pars[d1l:d1u,d2i,d3i]\
							= datap['m1par'][vl:vu]
							M2pars[d1l:d1u,d2i,d3i]\
							= datap['m2par'][vl:vu]
							M3pars[d1l:d1u,d2i,d3i]\
							= datap['m3par'][vl:vu]
						if(bGrav):
							Phis[d1l:d1u,d2i,d3i]\
							= datap['phi'][vl:vu]

				for ii in range(nxp):
					for kk in range(nzp): 
						d1i = ii+cnt1
						d2l = Nx2-cnt2-nyp
						d2u = Nx2-cnt2
						d3i = Nx3-cnt3-kk-1

						x2s[d1i,d2l:d2u,d3i] =\
						np.fliplr([datap['x2']])[0]			

				for ii in range(nxp):
					for jj in range(nyp):
						d1i = ii+cnt1
						d2i = Nx2-cnt2-jj-1
						d3l = Nx3-cnt3-nzp
						d3u = Nx3-cnt3

						x3s[d1i,d2i,d3l:d3u] =\
						np.fliplr([datap['x3']])[0]
				
				# advance positions in global 3D array indices
				cnt1 += nxp
			cnt2 += nyp
		cnt3 += nzp

	# flip 3D arrays for orientation in matplotlib imshow
	x1s = bin3Dflips(x1s)
	x2s = bin3Dflips(x2s)
	x3s = bin3Dflips(x3s)
	dps = bin3Dflips(dps)
	V1s = bin3Dflips(V1s)
	V2s = bin3Dflips(V2s)
	V3s = bin3Dflips(V3s)
	if(bPar):
		dpars  = bin3Dflips(dpars)
		M1pars = bin3Dflips(M1pars)
		M2pars = bin3Dflips(M2pars)
		M3pars = bin3Dflips(M3pars)
	if(bGrav):
		Phis = bin3Dflips(Phis)

	# Store global 3D arrays in a single dictionary, return that
	data = {}
	data['x1'] = x1s
	data['x2'] = x2s
	data['x3'] = x3s
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

	return data


def read_bin_header(filename,bDoublePres=True):
	"""Reads header of athena .bin file for time data.
	
	Can easily be movied to return other quantities if wanted.

	Parameters
	-------------
	filename : str
		full file name
	bDoublePres : bool, optional
		True if double precision was used in data output
	
	Returns
	-----------
	floats
		t, dt. Current sim. time and curr. time step.
	"""

	try:
	  file = open(filename,'rb')
	except:
	  print("Couldn't open file, tried: ", filename)
	  raise SystemExit
	print(filename)

	floatSize = np.float32
	filesize_float = 4
	if(bDoublePres == 1):
		floatSize = np.float64
		filesize_float = 8

	file.seek(0,2)
	eof = file.tell()
	file.seek(0,0)

	filesize = 0

	coordsys = np.fromfile(file,dtype=np.int32, count=1)[0]
	nx,ny,nz = np.fromfile(file,dtype=np.int32,count=3)[0:3]

	filesize += 4*4

	NVAR = np.fromfile(file,dtype=np.int32,count=1)[0]
	NSCALARS = np.fromfile(file,dtype=np.int32,count=1)[0]
	iSelfGravity = np.fromfile(file,dtype=np.int32,count=1)[0]
	iParticles = np.fromfile(file,dtype=np.int32,count=1)[0]

	filesize += 4*4

	gamma1,cs = np.fromfile(file,dtype=floatSize,count=2)
	t,dt = np.fromfile(file,dtype=floatSize,count=2)

	filesize += filesize_float * 4

	return t, dt

# -------------------------------------------------------------------------- #
# ------------------------- Internal functions ----------------------------- #
# -------------------------------------------------------------------------- #

def bin3Dflips(arr):
	"""Flip 3D arrays to orient with matplotlib axes.
	
	"""
	return np.flip(np.flip(arr,axis=1),axis=2)

def parse_bin(filename,bDoublePres):
	"""Reads data from athena .bin file and sets up data reading.
	
	"""

	try:
	  file = open(filename,'rb')
	except:
	  print("Couldn't open file, tried: \n", filename)
	  raise SystemExit
	print(filename)

	floatSize = np.float32
	filesize_float = 4
	if(bDoublePres):
		floatSize = np.float64
		filesize_float = 8

	file.seek(0,2)
	eof = file.tell()
	file.seek(0,0)

	filesize = 0

	coordsys = np.fromfile(file,dtype=np.int32, count=1)[0]
	nx,ny,nz = np.fromfile(file,dtype=np.int32,count=3)[0:3]

	filesize += 4*4

	NVAR = np.fromfile(file,dtype=np.int32,count=1)[0]
	NSCALARS = np.fromfile(file,dtype=np.int32,count=1)[0]
	iSelfGravity = np.fromfile(file,dtype=np.int32,count=1)[0]
	iParticles = np.fromfile(file,dtype=np.int32,count=1)[0]

	filesize += 4*4

	gamma1,cs = np.fromfile(file,dtype=floatSize,count=2)
	t,dt = np.fromfile(file,dtype=floatSize,count=2)

	filesize += filesize_float * 4

	x1 = np.fromfile(file,dtype=floatSize,count=nx)
	x2 = np.fromfile(file,dtype=floatSize,count=ny)
	x3 = np.fromfile(file,dtype=floatSize,count=nz)

	filesize += filesize_float * (nx+ny+nz)

	NVAR = NVAR.astype(int)
	nx = nx.astype(int)
	ny = ny.astype(int)
	nz = nz.astype(int)

	# Read regular variables
	data_read = np.zeros([NVAR,nz*ny*nx])
	for n in range(NVAR):
		data_read_3D, filesize = athRead3D(file, nz, ny, nx, floatSize, 
						filesize, filesize_float)

		data_read[n,:] = data_read_3D

	cell_d = data_read[0,:]
	cell_V1 = data_read[1,:]
	cell_V2 = data_read[2,:]
	cell_V3 = data_read[3,:]

	# If self gravity is on
	if( iSelfGravity == 1):
		data_read_3D, filesize = athRead3D(file, nz, ny, nx, floatSize, 						filesize, filesize_float)
		Phi = data_read_3D

	# If particles is on
	# Read particle variables
	if( iParticles == 1):

		data_read = np.zeros([4,nz*ny*nx]) # Hardcoded NVAR = 4,it seems
		for n in range(4):
			data_read_3D, filesize = athRead3D(file, nz, ny, nx,
					floatSize, filesize, filesize_float)
			data_read[n,:] = data_read_3D

		dpar = data_read[0,:]

		M1par = data_read[1,:] #/ dpar
		M2par = data_read[2,:] #/ dpar
		M3par = data_read[3,:] #/ dpar

	#print('Filesize is: ', filesize)

	data = {}
	data['x1'] = x1
	data['x2'] = x2
	data['x3'] = x3
	data['d'] = cell_d
	data['v1'] = cell_V1
	data['v2'] = cell_V2
	data['v3'] = cell_V3

	if(iParticles==1):
		data['dpar'] = dpar
		data['m1par'] = M1par
		data['m2par'] = M2par
		data['m3par'] = M3par
	if(iSelfGravity==1):
		data['phi'] = Phi

	return nx, ny, nz, data

def athRead3D(file, nz, ny, nx, floatSize, filesize, filesize_float):	
	"""Reads grid data from athena .bin file.
	
	"""

	data_read_3D = np.asarray([])
	data_read_once = np.zeros([nx])

	for k in range(nz):
		data_read_2D = np.asarray([])
	
		for j in range(ny):
				
			data_read_once = np.fromfile(file,dtype=floatSize,
							count=nx)

			filesize += filesize_float * nx


			data_read_2D = np.concatenate((data_read_2D,
							data_read_once))

		data_read_3D = np.concatenate((data_read_3D,data_read_2D))	

	return data_read_3D, filesize


