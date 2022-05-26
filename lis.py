import numpy as np
import struct

def read_lis(probid, filenum, filext, numprocs,
	bDoublePres=False):	
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
	bPar : bool, optional
		Boolean for if particle module was used
	bDoublePres : bool, optional
		True if double precision was used in data output

	Returns
	-----------
	dict
		3D numpy arrays
	"""

	px1 = np.asarray([])
	px2 = np.asarray([])
	px3 = np.asarray([])
	pv1 = np.asarray([])
	pv2 = np.asarray([])
	pv3 = np.asarray([])
	pdpar = np.asarray([])
	pdpar_St = np.asarray([])
	pgrproperty = np.asarray([])
	pmy_id = np.asarray([])
	pinit_id = np.asarray([])

	for idn3 in range(numprocs[2]):	
		cnt2 = 0;
		for idn2 in range(numprocs[1]):
			cnt1 = 0
			for idn1 in range(numprocs[0]):
				idno = idn1+numproc_x1*idn2+numproc_x1*\
					numproc_x2*idn3

				filepath = './data'
				afterPre = ''	
				if( numprocs[0] == 1 and
					numprocs[1] == 1 and
					numprocs[2] == 1):
					filepath += ''
					afterPre += ''
				else:
					if(idno == 0):
						filepath += '/id' + str(idno)
						afterPre += ''
					else:
						filepath += '/id' + str(idno)
						afterPre += '-id' + str(idno)
				filepath += '/'
				afterPre += '.'

				fileEnd = afterPre + filenum + filext
				filename = filepath + probid + fileEnd

				datan, grpropertyrad = \
					read_single_lis(filename, bDoublePres)

				px1n = datan['px1']
				px2n = datan['px2']
				px3n = datan['px3']
				pv1n = datan['pv1']
				pv2n = datan['pv2']
				pv3n = datan['pv3']
				pdparn = datan['pdpar']
				pgrpropertyn = datan['grpropertys']
				pmy_id = datan['my_ids']
				pinit_id = datan['init_ids']

				px1 = np.concatenate((px1,px1n))
				px2 = np.concatenate((px2,px2n))
				px3 = np.concatenate((px3,px3n))
				pv1 = np.concatenate((pv1,pv1n))
				pv2 = np.concatenate((pv2,pv2n))
				pv3 = np.concatenate((pv3,pv3n))
				pdpar = np.concatenate((pdpar,pdparn))
				pgrproperty = np.concatenate((pgrproperty,\
					pgrpropertyn))
				pmy_id = np.concatenate((pmy_id,pmy_idn))
				pinit_id = np.concatenate((pinit_id,pinit_idn))
		
	pmy_id = pmy_id.astype('int')
	pinit_id = pinit_id.astype('int')

	return px1, px2, px3, pv1, pv2, pv3, pdpar,\
		pgrproperty, pmy_id, pinit_id, grpropertyrad



def read_single_lis(filename, bDoublePres=False):
	"""Reads single .lis file.
	
	Useful if all .lis files have already been stitched together
	with join_lis.c routine provided in athena repo.

	Parameters
	-------------
	filename : str
		full file name
	bDoublePres : bool, optional
		True if double precision was used in data output
	
	Returns
	-----------
	dict
		dictionary of 3D numpy arrays.
	"""


	print(filename)

	try:
	  file = open(filename,'rb')
	except:
	  print("Couldn't open file.")
	  raise SystemExit

	floatSize = np.float32
	filesize_float = 4
	if(bDoublePres == 1):
		floatSize = np.float64
		filesize_float = 8

	file.seek(0,2)
	eof = file.tell()
	file.seek(0,0)

	filesize = 0

	GridData = np.fromfile(file,dtype=floatSize, count=12)

	#print(GridData)

	filesize += filesize_float*12

	npartypes = np.fromfile(file,dtype=np.int32, count=1)[0]
	filesize += 4
	#print('npartypes', npartypes)

	grpropertyrad = np.zeros(npartypes)
	for i in range(npartypes):
		grpropertyrad[i] = np.fromfile(file,dtype=floatSize, count=1)[0]
		filesize += filesize_float
		#print('grpropertyrad', grpropertyrad[i])

	t,dt = np.fromfile(file,dtype=floatSize,count=2)
	filesize += filesize_float*2
	#print('t, dt', t, dt)

	nout = np.fromfile(file,dtype=np.int64, count=1)[0]
	#print('nout', nout)
	filesize += 8

	parData  = np.zeros([nout,7])
	grpropertys   = np.zeros(nout)
	my_ids   = np.zeros(nout)
	init_ids = np.zeros(nout)

	# buffer size for single particle's data
	bufferSize = int(filesize_float*7 + 4 + 8 + 4)

	for ii in range(nout): # number of particles output to .lis file
		# unpack data for single particle
		px1u, px2u, px3u, pv1u, pv2u, pv3u, pdparu, \
			grpropertyu, my_idu, init_idu = \
			struct.unpack("<fffffffiqi", file.read(bufferSize))

		parData[ii,:] = [px1u, px2u, px3u, pv1u, pv2u, pv3u, pdparu]
		filesize += filesize_float*7

		grpropertys[ii] = grpropertyu
		filesize += 4

		my_ids[ii] = my_idu
		filesize += 8

		init_ids[ii] = init_idu
		filesize += 4

	my_ids   = my_ids.astype('int')
	init_ids = init_ids.astype('int')

	data = {}
	data['px1'] = parData[:,0]
	data['px2'] = parData[:,1]
	data['px3'] = parData[:,2]
	data['pv1'] = parData[:,3]
	data['pv2'] = parData[:,4]
	data['pv3'] = parData[:,5]

	# some estimate of the particle density, not recommended to be used
	data['pdpar'] = parData[:,6] 

	data['my_ids'] = my_ids # particle id no. on an mpi process
	data['init_ids'] = init_ids # mpi process id no
	data['grpropertys'] = grpropertys # integer for grain species

	return data, grpropertyrad



def read_lis_header(filename, bDoublePres=False):
	"""Reads header of athena .lis file for time data.
	
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
	  print("Couldn't open file.")
	  raise SystemExit

	floatSize = np.float32
	filesize_float = 4
	if(bDoublePres == 1):
		floatSize = np.float64
		filesize_float = 8

	file.seek(0,2)
	eof = file.tell()
	file.seek(0,0)

	filesize = 0

	GridData = np.fromfile(file,dtype=floatSize, count=12)

	#print(GridData)

	filesize += filesize_float*12

	npartypes = np.fromfile(file,dtype=np.int32, count=1)[0]
	filesize += 4
	#print('npartypes', npartypes)

	grpropertyrad = np.fromfile(file,dtype=floatSize, count=npartypes)[0]
	filesize += filesize_float
	#print('grpropertyrad', grpropertyrad)

	t,dt = np.fromfile(file,dtype=floatSize,count=2)

	return t, dt

