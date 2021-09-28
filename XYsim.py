import filemod as fm
from decimal import Decimal
import numpy as np

def modifyAthenaInput(filename, vec, val, Kxn, Kzn, mu, tau_s, omega, 
		etavkOncs, nwaveX, nwaveY, amp, tau0, tauf, Nout, Lx, Ly, Nx1, Nx2):

	try:
		f = open(filename,"a")
	except:
		print("Couldn't open file.")
		raise SystemExit

	f.write('Kxn          = %.12E' % Decimal(Kxn) + '\n')
	f.write('Kyn          = %.12E' % Decimal(Kzn) + '\n')

	f.write('Reux         = %.12E' % Decimal(np.real(vec[4])) + '\n')
	f.write('Imux         = %.12E' % Decimal(np.imag(vec[4])) + '\n')
	f.write('Reuy         = %.12E' % Decimal(np.real(vec[5])) +  '\n')
	f.write('Imuy         = %.12E' % Decimal(np.imag(vec[5])) +  '\n')
	f.write('Rerho        = %.12E' % Decimal(np.real(vec[3])) + '\n')
	f.write('Imrho        = %.12E' % Decimal(np.imag(vec[3])) + '\n')
	f.write('Rewx         = %.12E' % Decimal(np.real(vec[1])) + '\n')
	f.write('Imwx         = %.12E' % Decimal(np.imag(vec[1])) + '\n')
	f.write('Rewy         = %.12E' % Decimal(np.real(vec[2])) + '\n')
	f.write('Imwy         = %.12E' % Decimal(np.imag(vec[2])) + '\n')

	f.write('omg          = %.12E' % Decimal(np.real(val)) + '\n')
	f.write('s            = %.12E' % Decimal(np.imag(val)) + '\n')

	f.write('mratio       = %.12E' % Decimal(mu) + '\n')
	f.write('tstop        = %.12E' % Decimal(tau_s) + '\n')
	f.write('etavk        = %.12E' % Decimal(etavkOncs) + '\n')

	f.write('nwaveX       = %.12E' % Decimal(nwaveX) + '\n')
	f.write('nwaveY       = %.12E' % Decimal(nwaveY) + '\n')

	f.close()

	try:
		f = open(filename,"a")
	except:
		print("Couldn't open file.")
		raise SystemExit

	tlim = (2./(3.*omega))*(tauf-tau0)
	bindt = tlim/Nout

	fm.lineReplace(filename, 'tlim=1', 
					'tlim            = {:f}'.format(tlim))
	fm.lineReplace(filename, 'bindt=1', 
					'dt      = {:f}'.format(bindt))
	fm.lineReplace(filename, 'x1min=500', 
					'x1min           = {:f}'.format(-Lx/2.))
	fm.lineReplace(filename, 'x1max=500', 
					'x1max           =  {:f}'.format( Lx/2.))
	fm.lineReplace(filename, 'x2min=500', 
					'x2min           = {:f}'.format(-Ly/2.))
	fm.lineReplace(filename, 'x2max=500', 
					'x2max           =  {:f}'.format( Ly/2.))
	fm.lineReplace(filename, 'Nx1=500', 
					'Nx1             = {:d}'.format( Nx1  ))
	fm.lineReplace(filename, 'Nx2=500', 
					'Nx2             = {:d}'.format( Nx2  ))
	fm.lineReplace(filename, 'amp=500', 
					'amp             = %.12E' % Decimal(amp) )
	f.close()

