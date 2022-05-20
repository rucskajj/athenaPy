import numpy as np
import sys
from matplotlib import pyplot
import struct

def readAthenaBin(fileName, filepath, filenum, numprocs, Nx, bPar, bGrav):

	Nx1 = Nx[0]; Nx2 = Nx[1]; Nx3 = Nx[2];

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

				filePath = filepath
				afterPre = ''
		
				if( numproc_x1 == 1 and
					numproc_x2 == 1 and
					numproc_x3 == 1):
					filePath += ''
					afterPre += ''
				else:
					if(idno == 0):
						filePath += 'id' + str(idno)
						afterPre += ''
					else:
						filePath += 'id' + str(idno)
						afterPre += '-id' + str(idno)

				filePath += '/'
				afterPre += '.'

				fileEnd = afterPre + filenum + '.bin'

				filename = filePath + fileName + fileEnd

				print(filename)
	
				if(bPar==1 and bGrav==1):			
					nxp, nyp, nzp, x1p, x2p, x3p, dp, V1p, V2p, V3p, dparp,\
						M1parp,M2parp, M3parp, Phip = \
						parseAthenaBin(filename, bDoublePres)

				elif(bPar==1 and bGrav==0):			
					nxp, nyp, nzp, x1p, x2p, x3p, dp, V1p, V2p, V3p, dparp,\
						M1parp,M2parp, M3parp = \
						parseAthenaBin(filename, bDoublePres)

				elif(bPar==0 and bGrav==1):			
					nxp, nyp, nzp, x1p, x2p, x3p, dp, V1p, V2p, V3p, Phip = \
						parseAthenaBin(filename, bDoublePres)

				elif(bPar==0 and bGrav==0):			
					nxp, nyp, nzp, x1p, x2p, x3p, dp, V1p, V2p, V3p = \
						parseAthenaBin(filename, bDoublePres)


				#print('')
				#print(idno, ':', dparp)
				#print('')

				for ii in range(nxp):
					for kk in range(nzp): 
						x2s[ii+indCount1,(Nx2-indCount2-nyp):(Nx2-indCount2),Nx3-indCount3-kk-1] = np.fliplr([x2p])[0]			
					
				for kk in range(nzp):		
					for jj in range(nyp):
						x1s[indCount1:nxp+indCount1,Nx2-indCount2-jj-1,Nx3-indCount3-kk-1] = x1p
						dps[indCount1:nxp+indCount1,Nx2-indCount2-jj-1,Nx3-indCount3-kk-1] = dp[(jj*nxp)+(kk*nxp*nyp):(jj+1)*nxp+kk*nxp*nyp]
						V1s[indCount1:nxp+indCount1,Nx2-indCount2-jj-1,Nx3-indCount3-kk-1] = V1p[(jj*nxp)+(kk*nxp*nyp):(jj+1)*nxp+kk*nxp*nyp]
						V2s[indCount1:nxp+indCount1,Nx2-indCount2-jj-1,Nx3-indCount3-kk-1] = V2p[(jj*nxp)+(kk*nxp*nyp):(jj+1)*nxp+kk*nxp*nyp]
						V3s[indCount1:nxp+indCount1,Nx2-indCount2-jj-1,Nx3-indCount3-kk-1] = V3p[(jj*nxp)+(kk*nxp*nyp):(jj+1)*nxp+kk*nxp*nyp]

						if(bPar):
							dpars[indCount1:nxp+indCount1,Nx2-indCount2-jj-1,Nx3-indCount3-kk-1] = dparp[(jj*nxp)+(kk*nxp*nyp):(jj+1)*nxp+kk*nxp*nyp]
							M1pars[indCount1:nxp+indCount1,Nx2-indCount2-jj-1,Nx3-indCount3-kk-1] = M1parp[(jj*nxp)+(kk*nxp*nyp):(jj+1)*nxp+kk*nxp*nyp]
							M2pars[indCount1:nxp+indCount1,Nx2-indCount2-jj-1,Nx3-indCount3-kk-1] = M2parp[(jj*nxp)+(kk*nxp*nyp):(jj+1)*nxp+kk*nxp*nyp]
							M3pars[indCount1:nxp+indCount1,Nx2-indCount2-jj-1,Nx3-indCount3-kk-1] = M3parp[(jj*nxp)+(kk*nxp*nyp):(jj+1)*nxp+kk*nxp*nyp]


						if(bGrav):
							Phis[indCount1:nxp+indCount1,Nx2-indCount2-jj-1,Nx3-indCount3-kk-1] = Phip[(jj*nxp)+(kk*nxp*nyp):(jj+1)*nxp+kk*nxp*nyp]

						#print(dparp[(jj*nxp)+(kk*nxp*nyp):(jj+1)*nxp+kk*nxp*nyp])
						#print('')

				for ii in range(nxp):
					for jj in range(nyp):
						x3s[ii+indCount1,Nx2-indCount2-jj-1,(Nx3-indCount3-nzp):(Nx3-indCount3)] = np.fliplr([x3p])[0]



				indCount1 += nxp
			indCount2 += nyp
		indCount3 += nzp

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

	
	if(bPar==1 and bGrav==1):
		return x1s, x2s, x3s, dps, V1s, V2s, V3s, dpars, M1pars, M2pars, \
			M3pars, Phis

	elif(bPar==1 and bGrav==0):
		return x1s, x2s, x3s, dps, V1s, V2s, V3s, dpars, M1pars, M2pars, M3pars

	elif(bPar==0 and bGrav==1):
		return x1s, x2s, x3s, dps, V1s, V2s, V3s, Phis

	elif(bPar==0 and bGrav==0):
		return x1s, x2s, x3s, dps, V1s, V2s, V3s

def bin3Dflips(arr):
	return np.flip(np.flip(arr,axis=1),axis=2)

def transfer3Dto2D(A, Nx):
	Nx1 = Nx[0]; Nx2 = Nx[1]; Nx3 = Nx[2];

	B = np.zeros(Nx1*Nx2*Nx3)

	for kk in range(Nx3):
		for jj in range(Nx2):
			B[(jj*Nx1)+(kk*Nx2*Nx1):((jj+1)*Nx1)+(kk*Nx2*Nx1)] = A[:,jj,kk]

	return B

def parseAthenaBin(filename,bDoublePres):

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

	if(iParticles==1 and iSelfGravity==1):
		return nx, ny, nz, x1, x2, x3, cell_d, cell_V1, cell_V2, cell_V3, \
			dpar, M1par, M2par, M3par, Phi

	elif(iParticles==1 and iSelfGravity==0):
		return nx, ny, nz, x1, x2, x3, cell_d, cell_V1, cell_V2, cell_V3,\
			dpar, M1par, M2par, M3par

	elif(iParticles==0 and iSelfGravity==1):
		return nx, ny, nz, x1, x2, x3, cell_d, cell_V1, cell_V2, cell_V3, Phi

	elif(iParticles==0 and iSelfGravity==0):
		return nx, ny, nz, x1, x2, x3, cell_d, cell_V1, cell_V2, cell_V3


def athRead3D(file, nz, ny, nx, floatSize, filesize, filesize_float):	
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


def readAthenaBinHeader(filename,bDoublePres):

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

# ---------------------------------------------------------------------------- #

def readAthenaLis(fileName, filenum, filext, numprocs,bDoStitches):	
# for reading all *.lis files from the id*/ directories

	numproc_x1 = numprocs[0]
	numproc_x2 = numprocs[1]
	numproc_x3 = numprocs[2]

	bDoublePres = 0 

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

	for idn3 in range(numproc_x3):	
		indCount2 = 0;
		for idn2 in range(numproc_x2):
			indCount1 = 0
			for idn1 in range(numproc_x1):
				idno = idn1+numproc_x1*idn2+numproc_x1*\
					numproc_x2*idn3

				filePath = './data'
				afterPre = ''
		
				if( numproc_x1 == 1 and
					numproc_x2 == 1 and
					numproc_x3 == 1):
					filePath += ''
					afterPre += ''
				else:
					if(idno == 0):
						filePath += '/id' + str(idno)
						afterPre += ''
					else:
						filePath += '/id' + str(idno)
						afterPre += '-id' + str(idno)

				filePath += '/'
				afterPre += '.'

				fileEnd = afterPre + filenum + filext

				filename = filePath + fileName + fileEnd

				px1n, px2n, px3n, pv1n, pv2n, pv3n, pdparn, \
					pgrpropertyn, pmy_idn, pinit_idn,\
					t, dt, grpropertyrad = \
					readAthenaLisSingle(filename, \
						bDoublePres)

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



def readAthenaLisSingle(filename,bDoublePres=0):
# if join_lis has been used to put all .lis files together

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

	bufferSize = int(filesize_float*7 + 4 + 8 + 4)
	#print(bufferSize)

	for ii in range(nout): # number of particles output to .lis file
		px1u, px2u, px3u, pv1u, pv2u, pv3u, pdparu, \
			grpropertyL, my_idL, init_idL = \
			struct.unpack("<fffffffiqi", file.read(bufferSize))

		parData[ii,:] = [px1u, px2u, px3u, pv1u, pv2u, pv3u, pdparu]
		filesize += filesize_float*7

		grpropertys[ii] = grpropertyL
		filesize += 4

		my_ids[ii] = my_idL
		filesize += 8

		init_ids[ii] = init_idL
		filesize += 4

	px1 = parData[:,0]
	px2 = parData[:,1]
	px3 = parData[:,2]
	pv1 = parData[:,3]
	pv2 = parData[:,4]
	pv3 = parData[:,5]
	pdpar = parData[:,6]

	#px2 -= (GridData[2]+GridData[3])

	my_ids   = my_ids.astype('int')
	init_ids = init_ids.astype('int')

	return px1, px2, px3, pv1, pv2, pv3, pdpar, \
		grpropertys, my_ids, init_ids, t, dt, grpropertyrad



def readAthenaLisHeader(filename, bDoublePres=0):

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

