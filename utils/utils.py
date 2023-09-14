# utilities functions for run.py 

import math
import mesa_reader as mr
import numpy as np
import warnings
import os
import re
import sys

# check of whether scipy can be imported
try:
	from scipy.optimize import curve_fit
	from scipy.integrate import solve_ivp
	import scipy.interpolate as scipl
	scipy_exists = True
except:
	scipy_exists = False
	pass
# check of whether matplotlib can be imported
try:
	from matplotlib import pyplot as plt
	matplotlib_exists = True
except:
	matplotlib_exists = False
	pass

MSUN = 1.989E+33
RSUN = 6.963E+10
G = 6.6743e-8
yr = 3.15576E+07
pi = np.pi

# for CSM density profile fitting
def CSMprof_func(r, r_break, rho_break, yrho):
	nmax = 12.
	return np.log( rho_break * (( (np.exp(r) / r_break)**(1.5/yrho) + (np.exp(r) / r_break)**(nmax/yrho) ) /2. )**(-yrho) )

# for CSM density profile (Type Ibn/Icn)
def CSMprof_func2(r, r_break, rho_break, nout):
	yrho = 2.5
	return np.log( rho_break * (( (np.exp(r) / r_break)**(1.5/yrho) + (np.exp(r) / r_break)**(nout/yrho) ) /2. )**(-yrho) )

# find mesa model when mass eruption occurs
def find_mass_eruption(data_files, data_file_at_core_collapse, time_till_collapse):
	delta = np.zeros(len(data_files))
	for i, data_file in enumerate(data_files):
		delta[i] = abs(mr.MesaData(data_file_at_core_collapse).star_age - mr.MesaData(data_file).star_age - time_till_collapse)
	return data_files[np.argmin(delta)]


# obtain time from mass eruption till core-collapse
def get_mass_eruption_to_core_collapse(data_file_at_mass_eruption, data_file_at_core_collapse):
	data_me = mr.MesaData(data_file_at_mass_eruption)
	data_cc = mr.MesaData(data_file_at_core_collapse)
	# in years
	return data_cc.star_age - data_me.star_age


# remnant mass from fitting formulae of Schneider+20, arXiv:2008.08599
def remnant_from_CO(CO_core_mass):
	if CO_core_mass<6.357 or (CO_core_mass > 7.311 and CO_core_mass < 12.925):
		Mrem = 0.03357 * CO_core_mass + 1.31780
	elif (CO_core_mass > 6.357 and CO_core_mass < 7.311):
		Mrem = 10.**(-0.02466*CO_core_mass+1.28070)
	else:
		Mrem = 10.**(0.01940*CO_core_mass+0.98462)
	return Mrem


# extractor of envelope profile from the script 
def cc_param_extractor(data_file):
	data = mr.MesaData(data_file)
	r_star = data.photosphere_r * RSUN
	total_mass = data.star_mass
	CO_core_mass = data.c_core_mass
	return r_star, total_mass, CO_core_mass


# ejecta calculation script
def calculate_ejecta(data_file, file_at_cc, file_CSM):
	# extract total mass and CO core mass of the progenitor
	r_star, total_mass, CO_core_mass = cc_param_extractor(data_file)
	# extract density and pressure profile inside the star
	renv = np.loadtxt(file_at_cc, skiprows=1)[:,2]
	rhoenv = np.loadtxt(file_at_cc, skiprows=1)[:,4]
	penv =  np.loadtxt(file_at_cc, skiprows=1)[:,7]
	lgrhoenv = [math.log(rho) for i, rho in enumerate(rhoenv) if renv[i]<r_star]
	lgpenv = [math.log(p) for i, p in enumerate(penv) if renv[i]<r_star]
	# obtain polytropic index from fitting rho vs p at envelope. We fit with P = K*rho^(1+1/N), where N is the polytripic index.
	gamma = np.polyfit(lgrhoenv, lgpenv, 1)[0]
	Npol = 1./(gamma-1.)
	# Determine n using Matzner & McKee 1999 eq 25
	beta = 0.19
	n = (Npol+1.+3.*beta*Npol)/(beta*Npol)
	assert n>5, "Ejecta index should be >5, now %e" % n
	# set delta to 1
	delta = 1.0
	# obtain CSM mass from CSM file
	MrCSM = np.loadtxt(file_CSM, skiprows=1)[:,1]
	CSM_mass = (MrCSM[-1] - MrCSM[0]) / MSUN 
	remnant_mass = remnant_from_CO(CO_core_mass)
	Mej = total_mass - remnant_mass - CSM_mass
	assert Mej > 0.0
	print("Mej:%f Msun, n:%f, delta:%f, Mcsm:%f Msun" % (Mej, n, delta, CSM_mass), file=sys.stderr)
	return Mej, n, delta, CSM_mass


# remesh CSM for input to the light curve calculation
def remesh_CSM(rmax, CSM_in, CSM_out, data_file_at_mass_eruption, Ncell=1000, analytical_CSM=False, steady_wind='attach'):
	# copy first line
	with open(CSM_in, 'r') as fin:
		head = fin.readline().rstrip('\n')	
	# record values with edited mesh
	data = mr.MesaData(data_file_at_mass_eruption)
	CSM = np.loadtxt(CSM_in, skiprows=1)
	r_in = CSM[:,2]
	v_in = CSM[:,3]
	rho_in = CSM[:,4]
	rmin = max(r_in[0], 2.*data.photosphere_r*RSUN)
	v_edge = CSM[-1,3]
	X_edge = CSM[-1,5]
	Y_edge = CSM[-1,6]
	p_in = CSM[:,7]
	if steady_wind == 'RSGwind':
		if data.star_mdot > 0.0:
			# use the value from MESA
			wind_Mdot = -data.star_mdot * MSUN / 3.15e7
		else:
			# Nieuwenhuijzen & de Jager 90
			wind_Mdot = 5.6e-6 * (data.photosphere_L/1e5)**1.64 * (data.Teff/3500.)**(-1.61) * MSUN / 3.15e7 
		# wind velocity from galactic RSGs (Mauron+11, Appendix C)
		v_wind = 2e6 * (data.photosphere_L/1e5)**0.35
		wind_Mdot_vw = wind_Mdot / v_wind

	# find where the CSM contains artificial shock
	# We identify this to be the outermost radius where the pressure slope suddenly jumps to <-100.
	# this can be the radius where the CSM density becomes unreliable if there exists an artificial shock, or simply
	# can be the boundary of the star and the CSM.
	slope = [r/p_in[i+1]*(p_in[i+1]-p_in[i])/(r_in[i+1]-r_in[i]) for i,r in enumerate(r_in) if i<len(r_in)-1]
	# narrow down by also requiring outside of this to have a density profile of index ~-1.5.
	# get a "global density slope" since there exists local numerical fluctuations
	outward_global_slope = [r_in[i+10]/rho_in[i+10]*(rho_in[i+20]-rho_in[i+10])/(r_in[i+20]-r_in[i+10]) for i,r in enumerate(r_in) if i<len(r_in)-20]
	istop = max([i for i in range(len(slope[:-21])) if slope[i]<-100 and outward_global_slope[i+1]>-4.0 and outward_global_slope[i+1]<0.0])
	# one cell forward to not include the jumped cell
	istop += 1 
	rstop = r_in[istop]

	if analytical_CSM and scipy_exists:
		popt, pcov = curve_fit(CSMprof_func, np.log(r_in[istop:-2]), np.log(rho_in[istop:-2]), p0=[1e15,1e-15,2.0])
		(r_break, rho_break, yrho) = (popt[0], popt[1], popt[2])
	
	# generate output file of remeshed CSM
	r_out = np.logspace(math.log10(rmin*1.001), math.log10(rmax*1.001), Ncell)
	Mr_out = np.zeros(Ncell)
	v_out = np.zeros(Ncell)
	rho_out = np.zeros(Ncell)
	X_out = np.zeros(Ncell)
	Y_out = np.zeros(Ncell)
	last_Mr = 0.0
	Y_avrg = 0.0
	counter = 0
	for i, r in enumerate(r_out):
		if r < rstop:
			# we fix the density profile where the CSM density become unreliable due to artificial shocks.
			if analytical_CSM and scipy_exists:
				rho_out[i] = math.exp(CSMprof_func(math.log(r), r_break, rho_break, yrho))
			else:
				rho_out[i] = rho_in[istop] * (r/rstop)**(-1.5)
			# Mr is obtained later once rho and r are determined
			# for other ones use values at istop, since they are not important anyway
			v_out[i] = CSM[istop, 3]
			X_out[i] = CSM[istop, 5]
			Y_out[i] = CSM[istop, 6]
			last_Mr = CSM[istop, 1]
			counter += 1
		elif r < r_in[-1]:
			# initial Mr is that at istop.
			# NOTE this sets the zero point of the enclosed mass, which is presumably close
			# to the actual value. We don't really mind because we only use the difference of Mr
			# between 2 cells, but if we include the central star's gravity for calculation of
			# the light curve, this has to be revisited.
			# obtain rho by log interpolation, Mr and v by linear. Mr is not used anyway
			index = len([thisr for thisr in r_in if thisr < r])-1
			fraction = (r - r_in[index]) / (r_in[index+1] - r_in[index])
			v_out[i] = CSM[index, 3]*(1.-fraction) + CSM[index+1,3]*(fraction)
			rho_out[i] = CSM[index, 4]**(1.-fraction) * CSM[index+1,4]**(fraction) 
			Mr_out[i] = last_Mr + 4.*math.pi*r**2*rho_out[i]*(r-r_out[i-1])
			X_out[i] = CSM[index, 5]*(1.-fraction) + CSM[index+1,5]*(fraction)
			Y_out[i] = CSM[index, 6]*(1.-fraction) + CSM[index+1,6]*(fraction)
			# for next step
			last_Mr = Mr_out[i]
		else:
			# use a profile that connects to the wind profile with a Gaussian drop 
			if steady_wind == 'RSGwind':
				rho_out[i] = (rho_in[-1]-wind_Mdot_vw/(4.*math.pi*r_in[-1]**2)) * math.exp(1.-(r/r_in[-1])**2)  + wind_Mdot_vw / (4.*math.pi*r**2)
				v_out[i] = v_wind
			# or connect a wind profile to the edge of the erupted CSM profile
			elif steady_wind == 'attach':
				rho_out[i] = rho_in[-1] * (r_in[-1]/r)**2
				v_out[i] = v_edge 
			Mr_out[i] = last_Mr + 4.*math.pi*r**2*rho_out[i]*(r-r_out[i-1])
			# X, Y are the edge value
			X_out[i] = X_edge 
			Y_out[i] = Y_edge 
			# for next step
			last_Mr = Mr_out[i]
	
	# obtain Mr at r < rstop
	for i in range(counter, -1, -1):
		Mr_out[i] = Mr_out[i+1] - 4.*math.pi*r_out[i+1]**2*rho_out[i+1]*(r_out[i+1]-r_out[i])

	# save the CSM to an output file
	np.savetxt(CSM_out, np.transpose([list(range(1,Ncell+1)), Mr_out, r_out, v_out, rho_out, X_out, Y_out]), fmt=['%d','%.8e','%.8e','%.8e','%.8e','%.8e','%.8e'], header=head)
	# plot CSM profile
	if matplotlib_exists:
		plt.rcParams["font.size"] = 14
		plt.xscale('log')
		plt.yscale('log')
		plt.xlabel('radius [cm]')
		plt.ylabel('density [g cm$^{-3}$]')
		plt.xlim(0.5*r_out[0], r_out[-1])
		plt.plot(r_out, rho_out, label='remeshed')
		plt.plot(r_in, rho_in, linestyle='dashed', label='original')
		plt.legend(loc='upper right')
		plt.grid(linestyle='dotted')
		plt.tight_layout()
		plt.savefig('LCFiles/CSM_comparison.png')

	# extract Y_avrg, needed for opacity calculation
	Y_avrg = np.inner(Y_out[1:], np.diff(Mr_out))/ (Mr_out[-1]-Mr_out[0])
	return Y_avrg


# extract peak luminosity and rise time, defined as the time from (frac*L_peak) to L_peak
def extract_peak_and_rise_time(LC_file, frac):
	time = np.loadtxt(LC_file)[:,0]
	lum = np.loadtxt(LC_file)[:,1]
	peaktime = time[np.argmax(lum)]
	peakL = max(lum)
	if frac <  lum[0] / peakL:
		return ValueError("Parameter frac too small to obtain rise time.")
	risetime = peaktime - time[np.argmin([abs(L-peakL*frac) for i,L in enumerate(lum) if time[i]<peaktime])]
	if frac < lum[-1] / peakL:
		warnings.warn("Parameter frac too small to obtain the decay time accurately...")
	decaytime = time[np.argmax(lum) - 1 + np.argmin([abs(L-peakL*frac) for i,L in enumerate(lum) if time[i]>peaktime])] - peaktime
	print("peak: %e erg/s. rise time: %e days, decay time:%e days" % (peakL, risetime, decaytime), file=sys.stderr)


# calculate the lightcurve of mass eruption
def get_mass_eruption_lightcurve(outputFile):
        print('Calculate mass eruption lightcurve', file=sys.stderr)
        maxFileSize = 90
        time = np.zeros(maxFileSize)
        luminosity = np.zeros(maxFileSize)
        temperature = np.zeros(maxFileSize)

        # Get the abundance
        fileName = 'EruptionFiles/atCCSN.txt'
        if os.path.exists(fileName) == False:
                X = 0.7
                Y = 0.28
                Z = 1 - X - Y
        if os.path.exists(fileName) == True:
                temp = np.loadtxt(fileName, skiprows=1)
                X = temp[len(temp)-1,5]
                Y = temp[len(temp)-1,6]
                Z = 1 - X - Y
        print('X='+str(X)+' Y='+str(Y)+' Z='+str(Z), file=sys.stderr)
        print('time(day) luminosity(erg/s) temperature(K)', file=sys.stderr)

        for i in range(1,maxFileSize):
                # Get the hydrodynamical result
                fileName='EruptionFiles/result%s.txt' % str(i).zfill(2)
                if os.path.exists(fileName) == False:
                        break
                data = np.loadtxt(fileName, skiprows=2)

                # Get the time
                f = open(fileName, 'r')
                data2 = f.readline()
                f.close()
                regex = re.compile('[+-]?(?:\d+\.?\d*|\.\d+)(?:(?:[eE][+-]?\d+)|(?:\*10\^[+-]?\d+))?')
                match = regex.findall(data2)
                time[i] = match[1]

                # Calculate luminosity and effective temperautre
                opacity = np.zeros(len(data))
                for j in range(0,len(data)):
                        scatter = 0.20*(1+X)/((1+2.7e11*data[j][4]/(data[j][8]*data[j][8]))*(1+np.power(data[j][8]/4.5e8,0.86)))
                        molecular = 0.1*Z
                        negH = 1.1e-25*np.power(Z*data[j][4],0.5)*np.power(data[j][8],7.7)
                        kramers = 4e25*(1+X)*(Z+0.001)*data[j][4]*np.power(data[j][8],-3.5)
                        opacity[j] = molecular + 1/((1/negH)+(1/(kramers+scatter)))

                photosphere = 1
                for j in reversed(range(1,len(data))):
                        tau = 0
                        for k in range(j,len(data)):
                                tau = tau + data[k][4]*(data[k][1] - data[k-1][1])*opacity[k]
                        if tau < 0.666666667:
                                photosphere = j
                        if tau > 0.666666667:
                                break
                photosphere = photosphere - 1
                luminosity[i] = data[photosphere][9]*1e40
                temperature[i] = data[photosphere][8]
                if photosphere == 0:
                        luminosity[i] = 0
        with open(outputFile, mode = 'w') as f:
                for i in range (1,len(time)):
                        f.write('{:.5e}'.format(time[i]/86400)+' '+'{:.5e}'.format(luminosity[i])+' '+'{:.5e}'.format(temperature[i])+'\n')
        f.close()


# https://ui.adsabs.harvard.edu/abs/2012MNRAS.422...70H/abstract
def discriminantProgModel(progName):
	h = mr.MesaData(progName)
	mass = h.mass
	h1 = h.h1
	he4 = h.he4
	dm = abs(np.diff(mass))
	dm = np.append(dm, mass[-1])

	Mass_H = np.sum(dm*h1)
	Mass_He = np.sum(dm*he4)

	if Mass_H > 0.033:
		return 0
	elif Mass_He > 0.03:
		return 1
	else:
		return 2


def evolv_CSM(tinj):
	t_final = tinj*yr
	filename = './EruptionFiles/result99.txt'

	try:
		data = np.loadtxt(filename, skiprows=4)
	except:
		print('The calculation for eruption has not finished yet. Exit.')
		sys.exit()
	
	f = open(filename)
	t_0 = float((f.readline().split(' ')[4]))
	r_star = np.loadtxt('./EruptionFiles/InitForHydro.txt', skiprows=2, usecols=4)[-1]
	f.close()

	r = data[:,1]
	Mr= data[:,2]
	dm = data[:,3]
	rho = data[:,4]
	v = data[:,5]

	condition = (v**2/2-G*Mr/r > 0.)
	r = r[condition]
	Mr= Mr[condition]
	dm = dm[condition]
	rho = rho[condition]
	v = v[condition]
	l = len(r)

	t_f = np.empty(l)
	r_f = np.empty(l)
	v_f = np.empty(l)

	t_span = [t_0,t_final]

	for i in range(l):
		Mr0 = Mr[i]
		init = [r[i],v[i]]
		def func(t,x): 
			r_c,v_c = x
			t = t
			dr = v_c
			dv = -G*Mr0/r_c**2
			return [dr,dv]
    
		sol = solve_ivp(func,t_span,init,method='RK45',atol = 1e-8,rtol = 1e-10)

		t_f[i] = sol.t[-1]
		r_f[i] = sol.y[0,:][-1]
		v_f[i] = sol.y[1,:][-1]

	rho_f = dm/(4*pi/3)/np.gradient(r_f**3)

	flag = (t_f >= t_final) 

	rho_f = rho_f[flag]
	v_f = v_f[flag]
	r_f = r_f[flag]

	l = len(r_f)
	result = np.empty((l,3))

	arg = np.argsort(r_f)
	result[:,0] = r_f
	result[:,1] = v_f
	result[:,2] = rho_f

	return result


def remesh_evolv_CSM(tinj, rout, CSM_out, data_file_at_mass_eruption, Ncell=1000):
	pi = np.pi
	result = evolv_CSM(tinj)
	data = mr.MesaData(data_file_at_mass_eruption)
	rho = result[:,2]
	r = result[:,0]
	v = result[:,1]
	v_esc = np.sqrt(2.*G*data.star_mass*MSUN/data.photosphere_r/RSUN)
	Mr = np.empty(Ncell)
	spline_v = scipl.CubicSpline(r, v)
	nsize = len(r)
	(rmin, rmax) = (r[0], r[-1])
	nout_init = -(np.log(rho[nsize-10])-np.log(rho[nsize-20]))/(np.log(r[nsize-10])-np.log(r[nsize-20]))

	popt, pcov = curve_fit(CSMprof_func2, np.log(r), np.log(rho), p0 = [4e15, 1e-16, nout_init])
	(r_star, rho_star, nout) = (popt[0], popt[1], popt[2])
	r_remesh = np.logspace(np.log10(rmin*1.001), np.log10(rout*1.001), Ncell)
	rho_remesh = np.zeros(Ncell)
	v_remesh = np.zeros(Ncell)

	for i, x in enumerate(r_remesh):
		if x <= rmax:
			rho_remesh[i] = np.exp(CSMprof_func2(np.log(x), r_star, rho_star, nout))
			v_remesh[i] = spline_v(x)
		else:
			v_remesh[i] = v_esc
			rho_remesh[i] = rho_remesh[i-1]*(x/r_remesh[i-1])**(-2.0)
		if i == 0:
			Mr[i] = 0.0
		else:
			Mr[i] = Mr[i-1]+4.0*pi*r_remesh[i]**2.0*rho_remesh[i]*(r_remesh[i]-r_remesh[i-1])
	head = 'CSM profile created by scipy.curve_fit'
	np.savetxt(CSM_out, np.transpose([list(range(1,Ncell+1)), Mr, r_remesh, v_remesh, rho_remesh, data.h1[0]*np.ones(Ncell), data.he4[0]*np.ones(Ncell)]), fmt=['%d','%.8e','%.8e','%.8e','%.8e','%.8e','%.8e'], header=head)


def genAbundanceTable(data_file_at_mass_eruption):
	data = mr.MesaData(data_file_at_mass_eruption)
	h1 = data.h1[0]
	he4 = data.he4[0]
	c12 = data.c12[0]
	o16 = data.o16[0]
	with open('input/abundance/abundance_for_tablegen.txt', mode='w') as f:
		s = '{:.4e}\n{:.4e}\n{:.4e}\n{:.4e}\n'.format(h1,he4,c12,o16)
		f.write(s)
