import numpy as np
from math import exp, pi

bands = ['u', 'b', 'v', 'r', 'i']
h = 6.62607015E-27
c = 2.99792458E+10
k = 1.380649E-16
a = 7.57E-15
pc = 3.085677581E+18

def B(nu, T):
	fac1 = 2.0*h/c**2.0*nu**3.0
	fac2 = h*nu/(k*T)
	if fac2 < 1.0e-05:
		return fac1/(fac2+fac2*(0.5*fac2+(fac2**2.0)/6.0))
	else:
		return fac1/(exp(fac2)-1.0)


	

def calc_mag(path):
	data = np.loadtxt(path)
	n = len(data[:,0])
	b = np.empty((n,len(bands)+1))
	b[:,0] = data[:,0]
	Tc = data[:,3]
	Lbol = data[:,1]
	r = np.sqrt(Lbol/(pi*a*c*Tc**4.0))

	for i, band in enumerate(bands):
		filename = './input/band_filters/'+band+'band.txt'
		fltr = np.loadtxt(filename)
		fltr[:,0] = fltr[:,0]*1.0e-07
		sum1 = np.zeros(n)
		sum2 = np.zeros(n)
		for j in range(len(fltr[:,0])-1):
			dlam = fltr[j+1,0]-fltr[j,0]
			lamm = 0.5*(fltr[j,0]+fltr[j+1,0])
			Trans = 0.5*(fltr[j,4]+fltr[j+1,4])
			array1 = [(r[l]/(10.0*pc))**2.0*pi*(B(c/fltr[j,0], Tc[l])+B(c/fltr[j+1,0], Tc[l]))/2.0*Trans*dlam/lamm for l in range(n)]
			array2 = [Trans*dlam/lamm for l in range(n)]
			sum1 += array1
			sum2 += array2
		
		b[:,i+1] = -2.5*np.log10(sum1/sum2)-48.6
	
	path = path.replace('.txt', '_mag.txt')
	np.savetxt(path, b)


###########################################################
#メモ帳：
#L_nu = A*B_nu(T_color)
#int_ L_nu dnu = L_bol
#int_ L_nu dnu = Asigma / pi * T_color^4
#でやっておく。
