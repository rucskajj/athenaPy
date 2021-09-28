import numpy as np
import sys
from matplotlib import pyplot

def readPeakTXT(filename):
	f = open(filename, 'r')

	header = f.readline()

	peak_id = []; Npar = [];
	total_mass = []; Hill_radius = [];
	com_x = []; com_y = []; com_z = [];
	velcom_x = []; velcom_y = []; velcom_z = [];
	Jx = []; Jy = []; Jz = [];

	for line in f:
		line = line.strip()
		columns = line.split()
		peak_id.append(int(columns[0]))
		Npar.append(int(columns[1]))
		total_mass.append(float(columns[2]))
		Hill_radius.append(float(columns[3]))
		com_x.append(float(columns[4]))
		com_y.append(float(columns[5]))
		com_z.append(float(columns[6]))
		velcom_x.append(float(columns[7]))
		velcom_y.append(float(columns[8]))
		velcom_z.append(float(columns[9]))
		Jx.append(float(columns[10]))
		Jy.append(float(columns[11]))
		Jz.append(float(columns[12]))

	peak_id = np.asarray(peak_id)
	Npar = np.asarray(Npar)
	total_mass = np.asarray(total_mass)
	Hill_radius = np.asarray(Hill_radius)
	com_x = np.asarray(com_x)
	com_y = np.asarray(com_y)
	com_z = np.asarray(com_z)
	velcom_x = np.asarray(velcom_x)
	velcom_y = np.asarray(velcom_y)
	velcom_z = np.asarray(velcom_z)
	Jx = np.asarray(Jx)
	Jy = np.asarray(Jy)
	Jz = np.asarray(Jz)

	return peak_id, Npar, total_mass, Hill_radius, com_x, com_y, com_z,\
		velcom_x, velcom_y, velcom_z, Jx, Jy, Jz

def readdparTXT(filename):
	f = open(filename, 'r')

	parid = []; 
	nden = []; aden = []; mass = [];
	x = []; y = []; z = [];
	velx = []; vely = []; velz = [];


	for line in f:
		line = line.strip()
		columns = line.split()
		parid.append(int(columns[0]))
		nden.append(float(columns[1]))
		aden.append(float(columns[2]))

		mass.append(float(columns[3]))

		x.append(float(columns[4]))
		y.append(float(columns[5]))
		z.append(float(columns[6]))

		velx.append(float(columns[7]))
		vely.append(float(columns[8]))
		velz.append(float(columns[9]))

	parid = np.asarray(parid)
	nden = np.asarray(nden)
	aden = np.asarray(aden)

	mass = np.asarray(mass)

	x = np.asarray(x)
	y = np.asarray(y)
	z = np.asarray(z)

	velx = np.asarray(velx)
	vely = np.asarray(vely)
	velz = np.asarray(velz)


	return parid, nden, aden, mass, x, y, z, velx, vely, velz

def readNaiveTXT(filename):
	f = open(filename, 'r')

	header = f.readline()

	x = []; y = []; z = [];
	dis_max = []; Npar = [];
	R_1_10 = []; R_halfM = []; R_moreM = [];


	for line in f:
		line = line.strip()
		columns = line.split()
		x.append(float(columns[0]))
		y.append(float(columns[1]))
		z.append(float(columns[2]))

		dis_max.append(float(columns[3]))
		Npar.append(int(columns[4]))

		R_1_10.append(float(columns[5]))
		R_halfM.append(float(columns[6]))
		R_moreM.append(float(columns[7]))

	x = np.asarray(x)
	y = np.asarray(y)
	z = np.asarray(z)

	dis_max = np.asarray(dis_max)
	Npar = np.asarray(Npar)

	R_1_10 = np.asarray(R_1_10)
	R_halfM = np.asarray(R_halfM)
	R_moreM = np.asarray(R_moreM)

	return x, y, z, dis_max, Npar, R_1_10, R_halfM, R_moreM
		
