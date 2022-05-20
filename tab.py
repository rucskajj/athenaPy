import numpy as np
import sys

# ---------------------------------------------------------------------------- #
# ----------------------------- 3D routines ---------------------------------- #
# ---------------------------------------------------------------------------- #

def readAthenaTab_3D(filepref, filePath, filenum, dataname, numprocs, 
	Nx, bDim):

	Nx1 = Nx[0]; Nx2 = Nx[1]; Nx3 = Nx[2];

	i = np.zeros((Nx1,Nx2,Nx3))
	j = np.zeros((Nx1,Nx2,Nx3))
	k = np.zeros((Nx1,Nx2,Nx3))
	datas = np.zeros((Nx1,Nx2,Nx3))
	
	numproc_x1 = numprocs[0]
	numproc_x2 = numprocs[1]
	numproc_x3 = numprocs[2]
	bDoublePres = 1

	indCount3 = 0
	for idn3 in range(numproc_x3):	
		indCount2 = 0
		for idn2 in range(numproc_x2):
			indCount1 = 0
			for idn1 in range(numproc_x1):

				idno = idn1+numproc_x1*idn2+numproc_x1*numproc_x2*idn3

				filepath = filePath
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

				fileEnd = afterPre + filenum + '.' + dataname + '.tab'
				fileEndbin = afterPre + '0000' + '.bin' 
			# Won't necessarily have a bin at every time step of the output, will have bin at 0000 always, I think.

				filename = filepath + filepref + fileEnd
				filenamebin = filepath + filepref + fileEndbin

				print(filename)

				ns = parseBin_for_ns(filenamebin, bDoublePres)
				nxp = ns[0]; nyp = ns[1]; nzp = ns[2];
				
				ip, jp, kp, datap = parseAthenaUsrExpr3D(filename, ns)

				for kk in range(nzp):		
					for jj in range(nyp):
						datas[indCount1:nxp+indCount1,Nx2-indCount2-jj-1,Nx3-indCount3-kk-1] = datap[(jj*nxp)+(kk*nxp*nyp):(jj+1)*nxp+kk*nxp*nyp]

				indCount1 += nxp
			indCount2 += nyp
		indCount3 += nzp
	
	return datas


def parseAthenaUsrExpr3D(filename, ns):

	try:
	  file = open(filename,'rb')
	except:
	  print("Couldn't open file.")
	  raise SystemExit

	i = np.zeros(ns[0]*ns[1]*ns[2])
	j = np.zeros(ns[0]*ns[1]*ns[2])
	k = np.zeros(ns[0]*ns[1]*ns[2])
	data = np.zeros(ns[0]*ns[1]*ns[2])

	ii = 0
	for line in file:
		dataline = np.asarray((line.strip()).split()).astype('float64')
		i[ii] = dataline[0].astype('int')
		j[ii] = dataline[1].astype('int')
		k[ii] = dataline[2].astype('int')
		data[ii] = dataline[3]
		ii += 1 

	return i, j, k, data

# ---------------------------------------------------------------------------- #
# ----------------------------- 2D routines ---------------------------------- #
# ---------------------------------------------------------------------------- #

def readTabs2D(probid, datapath, filenum, outid, Nx, xs, iDim, 
	numprocs=[1,1,1], bSlice=False):
	"""Reads 2D slice projection outputs from athena simulations.
	
	Parameters
	----------------

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
	xs : array_like
		output id from athena <output> block input block
	iDim: int
		what dimesion the slice or projection is along, values
		of 1, 2, 3 to match athena's x1, x2, x3 designation
	numprocs : array_like
		output id from athena <output> block input block
	bSlice : bool
		output id from athena <output> block input block

	
	Returns
	-------
	array_like
		2D slice of simulation data

	"""


	Nx1 = Nx[0]; Nx2 = Nx[1]; Nx3 = Nx[2];

	i = np.zeros((Nx1,Nx2,Nx3))
	j = np.zeros((Nx1,Nx2,Nx3))
	k = np.zeros((Nx1,Nx2,Nx3))
	datas = np.zeros((Nx1,Nx2,Nx3))
	
	numproc_x1 = numprocs[0]
	numproc_x2 = numprocs[1]
	numproc_x3 = numprocs[2]
	bDoublePres = 1

	filelist = []

	bFound = False
	indCount3 = 0
	for idn3 in range(numproc_x3):	
		indCount2 = 0
		for idn2 in range(numproc_x2):
			indCount1 = 0
			for idn1 in range(numproc_x1):

				idno = idn1+numproc_x1*idn2+numproc_x1*numproc_x2*idn3

				filepath = datapath
				# need to reset filepath to string from input
				# arg. list

				probiddir = ''
		
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
						filepath += 'id' + str(idno) + '/'
						probiddir += '-id' + str(idno) + '.'

				fileEnd = probiddir + filenum + '.' + outid + '.tab'
				filename = filepath + probid + fileEnd

				'''
				# Use .0000.bin files to fetch number of grid points in each
				# tab file per processor
				fileEndbin = probiddir + '0000' + '.bin' 
				filenamebin = filepath + probid + fileEndbin
				ns = parseBin_for_ns(filenamebin, bDoublePres)
				nxp = ns[0]; nyp = ns[1]; nzp = ns[2];
				'''

				try:
					file = open(filename,'rb')
					#print(filename)
				except:
					continue

				#parsedata = parseSingleTab2D(filename, ns)
				#print(parsedata)
				#if( parsedata is None ):
				#	continue

				filelist.append(filename)
				bFound = True
				file.close()

				#print(int(len(np.unique(parsedata[0]))), 
				#	int(len(np.unique(parsedata[1]))),ns)

				#n1p = int(len(np.unique(ip)))
				#n2p = int(len(np.unique(jp)))
				#print(n1p, n2p, len(datap), n1p*n2p)				

				'''

				for kk in range(nzp):		
					for jj in range(nyp):
						datas[indCount1:nxp+indCount1,Nx2-indCount2-jj-1,Nx3-indCount3-kk-1] = datap[(jj*nxp)+(kk*nxp*nyp):(jj+1)*nxp+kk*nxp*nyp]
				'''

#				indCount1 += nxp
#			indCount2 += nyp
#		indCount3 += nzp


	if(bFound == False):
		print("Couldn't open file, last tried: \n", filename)
		raise SystemExit	


	#print(Nx)
	inds = [1,2,3]
	inds.remove(iDim)
	inds[:] = [i-1 for i in inds]
	outdata = np.zeros([Nx[inds[0]], Nx[inds[1]]])
	print('outdata shape: ', outdata.shape)

	for filename in filelist:
		print(filename.replace(datapath, ''))

		binfilename = filename.replace(outid + '.tab', 'bin')
		binfilename = binfilename.replace(filenum, '0000')

		ns, ixs = parseBin_for_ns(binfilename, bDoublePres, xs)

		# -1 in min calculations for indexing numpy array later
		minx = ixs[inds[0]].min()-1; maxx = ixs[inds[0]].max();
		miny = ixs[inds[1]].min()-1; maxy = ixs[inds[1]].max();

		#print(minx, maxx, miny, maxy)
		#print(maxx-minx, maxy-miny)
		outdatabuf = np.zeros([maxx-minx, maxy-miny])

		#print(outdatabuf)
		#print('')
		outdatabuf = parseSingleTab2D(filename, outdatabuf)
		#print(outdatabuf)

		#print(minx, maxx, miny, maxy)
		#testchunk = outdata[minx:maxx, miny:maxy]
		#print(testchunk.shape, outdatabuf.shape)

		outdata[minx:maxx, miny:maxy] = outdatabuf

	return outdata



def parseSingleTab2D(filename, d):

	file = open(filename,'rb')

	ii = []
	jj = []
	data = []

	for line in file:
		dataline = np.asarray((line.strip()).split()).astype('float64')

		i = dataline[0].astype('int')
		j = dataline[1].astype('int')

		d[i,j] = dataline[2]

		ii.append(dataline[0].astype('int'))
		jj.append(dataline[1].astype('int'))
		data.append(dataline[2])


		# Using .append() method here as we may not know the number of grid
		# points given to each processor a priori

	file.close()

	return d


def parseBin_for_ns(filename, bDoublePres, xs):

	try:
		file = open(filename,'rb')
	except:
		print("Couldn't find 0000.bin file, tried: \n", filename)
		raise SystemExit

	floatSize = np.float32
	filesize_float = 4
	if(bDoublePres == 1):
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
	ix1 = np.round((x1+0.5*dx1 - xs[0][0])/dx1).astype('int')
	# ix1 is now global/domain wide cell index

	dx2 = x2[1]-x2[0]
	ix2 = np.round((x2+0.5*dx2 - xs[0][1])/dx2).astype('int')
	dx3 = x3[1]-x3[0]
	ix3 = np.round((x3+0.5*dx3 - xs[0][2])/dx3).astype('int')

	file.close()

	return [nx, ny, nz], [ix1, ix2, ix3]





