import numpy as np
import sys
from matplotlib import pyplot

def readAthenaUsrExprTab(filepref, filePath, filenum, dataname, numprocs, 
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

				if(bDim==1):
					nzp = 1 # If integrated along z

				if(bDim==0): #3D
					ip, jp, kp, datap = parseAthenaUsrExpr3D(filename, ns)
				elif(bDim==1): #2D
					ip, jp, datap = parseAthenaUsrExpr2D(filename, ns)

				#print('')
				#print(idno, ':', jp)
				#print('')

				#print(datap)

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

def parseAthenaUsrExpr2D(filename, ns):

	try:
	  file = open(filename,'rb')
	except:
	  print("Couldn't open file.")
	  raise SystemExit

	i = np.zeros(ns[0]*ns[1])
	j = np.zeros(ns[0]*ns[1])
	data = np.zeros(ns[0]*ns[1])

	ii = 0
	for line in file:
		dataline = np.asarray((line.strip()).split()).astype('float64')
		i[ii] = dataline[0].astype('int')
		j[ii] = dataline[1].astype('int')
		data[ii] = dataline[2]
		ii += 1 

	return i, j, data


def parseBin_for_ns(filename, bDoublePres):

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

	nx = nx.astype(int)
	ny = ny.astype(int)
	nz = nz.astype(int)

	return nx, ny, nz
