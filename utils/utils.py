# utilities functions for run.py 

import math
import mesa_reader as mr
import numpy as np
import warnings
import os
import re
import sys

from scipy.optimize import curve_fit
from scipy.integrate import solve_ivp, simpson
from scipy.interpolate import griddata
import scipy.interpolate as scipl

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
yr_to_sec = 3.15576E+07
pi = np.pi

# for CSM density profile fitting
def CSMprof_func(r, r_break, rho_break, yrho):
	nmax = 12.
	return np.log( rho_break * (( (np.exp(r) / r_break)**(1.5/yrho) + (np.exp(r) / r_break)**(nmax/yrho) ) /2. )**(-yrho) )

# for arbitrary CSM density profile fitting
def CSMprof_func_arb(r, r_break, rho_break, nin, yrho):
	nmax = 12.
	return rho_break * ((( r / r_break)**(-nin/yrho) + ( r / r_break)**(nmax/yrho) ) /2. )**(-yrho) 

# for CSM density profile (Type Ibn/Icn)
def CSMprof_func_stripped(r, r_break, rho_break, nout):
	yrho = 2.5
	return np.log( rho_break * (( (np.exp(r) / r_break)**(1.5/yrho) + (np.exp(r) / r_break)**(nout/yrho) ) /2. )**(-yrho) )

# find mesa model when mass eruption occurs
def find_mass_eruption(data_files, data_file_at_core_collapse, time_till_collapse):
	delta = np.zeros(len(data_files))
	for i, data_file in enumerate(data_files):
		delta[i] = abs(mr.MesaData(data_file_at_core_collapse).star_age - mr.MesaData(data_file).star_age - time_till_collapse)
	return data_files[np.argmin(delta)]


# remnant mass from fitting formulae of Schneider+20, arXiv:2008.08599 (single stars)
def remnant_from_CO_single_stars(CO_core_mass):
	if CO_core_mass<6.357 or (CO_core_mass > 7.311 and CO_core_mass < 12.925):
		Mrem = 0.03357 * CO_core_mass + 1.31780
	elif (CO_core_mass >= 6.357 and CO_core_mass < 7.311):
		Mrem = 10.**(-0.02466*CO_core_mass+1.28070)
	else:
		Mrem = 10.**(0.01940*CO_core_mass+0.98462)
	return Mrem

# remnant mass from fitting formulae of Schneider+20 for case B
def remnant_from_CO_case_B(CO_core_mass):
	if CO_core_mass<7.548 or (CO_core_mass > 8.491 and CO_core_mass < 15.144):
		Mrem = 0.01909 * CO_core_mass + 1.34529
	elif (CO_core_mass >= 7.548 and CO_core_mass < 8.491):
		Mrem = 10.**(0.03306*CO_core_mass+0.68978)
	else:
		Mrem = 10.**(0.02477*CO_core_mass+0.80614)
	return Mrem

# extractor of envelope profile from the script 
def cc_param_extractor(data_file, discriminant):
	data = mr.MesaData(data_file)
	r_star = data.photosphere_r * RSUN
	total_mass = data.star_mass
	# "c_core_mass" has been updated to "co_core_mass" from MESA r21.12.1
	try:
		CO_core_mass = data.c_core_mass
	except:
		CO_core_mass = data.co_core_mass
	# envelope mass coordinate where we measure the ejecta profile
	# We regard "envelope" as from 0.2Msun outside the "core" (in MESA), as how we do in the eruption calculations (see convert.py)
	if discriminant == 0:
		mr_env = data.he_core_mass+0.2
	elif discriminant == 1:
		mr_env = CO_core_mass+0.2
	elif discriminant == 2:
		mr_env = data.si_core_mass+0.2
	return r_star, total_mass, CO_core_mass, mr_env


# ejecta calculation script
def calculate_ejecta(data_file, file_at_cc, file_CSM, D):
	# extract total mass and CO core mass of the progenitor
	r_star, total_mass, CO_core_mass, mr_env = cc_param_extractor(data_file, D)
	# obtain ejecta outer power-law index n
	if 'atCCSN.txt' in file_at_cc:
		# obtaining n from simulated envelope at core-collapse
		renv = np.loadtxt(file_at_cc, skiprows=1)[:,2]
		rhoenv = np.loadtxt(file_at_cc, skiprows=1)[:,4]
		penv =  np.loadtxt(file_at_cc, skiprows=1)[:,7]
		lgrhoenv = [math.log(rho) for i, rho in enumerate(rhoenv) if renv[i]<r_star]
		lgpenv = [math.log(p) for i, p in enumerate(penv) if renv[i]<r_star]
	else:
		# obtaining n from envelope from MESA model
		h = mr.MesaData(data_file)
		lgrhoenv = h.logRho[h.mass > mr_env]
		lgpenv = h.logP[h.mass > mr_env]
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
	if D == 0:
		remnant_mass = remnant_from_CO_single_stars(CO_core_mass)
	elif D == 1 or D == 2:
		# assuming that the helium star losts its envelope by Case B mass transfer (which is most likely)
		remnant_mass = remnant_from_CO_case_B(CO_core_mass)
		
	Mej = total_mass - remnant_mass - CSM_mass
	assert Mej > 0.0
	print("Mej:%f Msun, n:%f, delta:%f, Mcsm:%f Msun" % (Mej, n, delta, CSM_mass), file=sys.stderr)
	return Mej, n, delta, CSM_mass


# remesh CSM for input to the light curve calculation
def remesh_CSM(rmax, CSM_in, CSM_out, data_file_at_mass_eruption, Ncell=1000, analytical_CSM=False, steady_wind='RSGwind'):
	# copy first line
	data = mr.MesaData(data_file_at_mass_eruption)
	# obtain (low-density) wind parameters, if we set this to be outside the dense CSM (default)
	if steady_wind == 'RSGwind':
		if data.star_mdot > 0.0:
			# use the value from MESA
			wind_Mdot = -data.star_mdot * MSUN / yr_to_sec
		else:
			# Nieuwenhuijzen & de Jager 90
			wind_Mdot = 5.6e-6 * (data.photosphere_L/1e5)**1.64 * (data.Teff/3500.)**(-1.61) * MSUN / yr_to_sec
		# wind velocity from galactic RSGs (Mauron+11, Appendix C)
		v_wind = 2e6 * (data.photosphere_L/1e5)**0.35
		wind_Mdot_vw = wind_Mdot / v_wind

	# if CSM from eruptive mass loss is used...
	if type(CSM_in) is str:
		with open(CSM_in, 'r') as fin:
			head = fin.readline().rstrip('\n')
		# record values with edited mesh
		CSM = np.loadtxt(CSM_in, skiprows=1)
		r_in = CSM[:,2]
		v_in = CSM[:,3]
		rho_in = CSM[:,4]
		rmin = max(r_in[0], 2.*data.photosphere_r*RSUN)
		v_edge = CSM[-1,3]
		X_edge = CSM[-1,5]
		Y_edge = CSM[-1,6]
		p_in = CSM[:,7]
	
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
	
		if analytical_CSM:
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
				if analytical_CSM:
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
				Mr_out[i] = last_Mr + 4.*pi*r**2*rho_out[i]*(r-r_out[i-1])
				X_out[i] = CSM[index, 5]*(1.-fraction) + CSM[index+1,5]*(fraction)
				Y_out[i] = CSM[index, 6]*(1.-fraction) + CSM[index+1,6]*(fraction)
				# for next step
				last_Mr = Mr_out[i]
			else:
				# use a profile that connects to the wind profile with a Gaussian drop 
				if steady_wind == 'RSGwind':
					rho_out[i] = (rho_in[-1]-wind_Mdot_vw/(4.*pi*r_in[-1]**2)) * math.exp(1.-(r/r_in[-1])**2)  + wind_Mdot_vw / (4.*pi*r**2)
					v_out[i] = v_wind
				# or connect a wind profile to the edge of the erupted CSM profile
				elif steady_wind == 'attach':
					rho_out[i] = rho_in[-1] * (r_in[-1]/r)**2
					v_out[i] = v_edge 
				Mr_out[i] = last_Mr + 4.*pi*r**2*rho_out[i]*(r-r_out[i-1])
				# X, Y are the edge value
				X_out[i] = X_edge 
				Y_out[i] = Y_edge 
				# for next step
				last_Mr = Mr_out[i]
		
		# obtain Mr at r < rstop
		for i in range(counter, -1, -1):
			Mr_out[i] = Mr_out[i+1] - 4.*pi*r_out[i+1]**2*rho_out[i+1]*(r_out[i+1]-r_out[i])
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

	else:
		print('The density profile of CSM is given by hand.')
		s = CSM_in[0]
		M_CSM = CSM_in[1]*MSUN
		r_break = CSM_in[2]
		rmin = data.photosphere_r*RSUN*2.
		if s <= -3.:
			print('power-law index should be smaller than 3.')
			sys.exit(1)
		print(s, M_CSM, r_break)

		r_out = np.logspace(math.log10(rmin*1.001), math.log10(rmax*1.001), Ncell)

		yrho = 3.
		rho_break = 1.e-10
		rho_break = M_CSM/(simpson(4.*pi*r_out**2.*CSMprof_func_arb(r_out, r_break, rho_break, s, yrho), x=r_out)/rho_break)
		rho_out = CSMprof_func_arb(r_out, r_break, rho_break, s, yrho)
		# connect RSG wind if requested
		if steady_wind == 'RSGwind':
			rho_out += wind_Mdot_vw / (4.*pi*r_out**2)
		# assume velocity is zero, and abundance is the surface abundance
		v_out = np.zeros(Ncell)
		X_out = np.ones(Ncell)*data.h1[0]
		Y_out = np.ones(Ncell)*data.he4[0]
		Mr_out = np.zeros(Ncell)
		dMr = 4.*pi*r_out**2*rho_out*np.gradient(r_out)
		for i in range(1, Ncell):
			Mr_out[i] = Mr_out[i-1]+dMr[i]
		print(Mr_out)
		head = 'The density profile of CSM created using arguments given by user'


	# save the CSM to an output file
	np.savetxt(CSM_out, np.transpose([list(range(1,Ncell+1)), Mr_out, r_out, v_out, rho_out, X_out, Y_out]), fmt=['%d','%.8e','%.8e','%.8e','%.8e','%.8e','%.8e'], header=head)

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


def discriminantProgModel(progName):
	h = mr.MesaData(progName)
	mass = h.mass
	h1 = h.h1
	he4 = h.he4
	dm = abs(np.diff(mass))
	dm = np.append(dm, mass[-1])

	# outer layer H mass
	Mass_H = np.sum(dm[mass>h.he_core_mass]*h1[mass>h.he_core_mass])

	# II/Ib border: https://ui.adsabs.harvard.edu/abs/2012MNRAS.422...70H/abstract
	if Mass_H > 0.033:
		return 0, 'IIn'
	else:
		# outer layer He mass
		Mass_He = np.sum(dm[mass>h.si_core_mass]*he4[mass>h.si_core_mass])
		# Ib/Ic border: https://ui.adsabs.harvard.edu/abs/2021ApJ...908..150W/abstract
		if Mass_He > 0.05:
			return 1, 'Ibn'
		else:
			return 2, 'Icn'


def evolv_CSM(tinj):
	t_final = tinj * yr_to_sec
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

	result[:,0] = r_f
	result[:,1] = v_f
	result[:,2] = rho_f

	return result


def remesh_evolv_CSM(tinj, rout, CSM_out, data_file_at_mass_eruption, Ncell=1000):
	result = evolv_CSM(tinj)
	data = mr.MesaData(data_file_at_mass_eruption)
	rho = result[:,2]
	r = result[:,0]
	v = result[:,1]
	v_esc = np.sqrt(2.*G*data.star_mass*MSUN/data.photosphere_r/RSUN)
	Mr = np.empty(Ncell)
	try:
		spline_v = scipl.CubicSpline(r, v)
	except: # there is a radius inversion in some of the cells.
		# remove that region, and let the interpolator take care of interpolation with the surrounding sane cells
		warnings.warn("Radius inversion in the CSM cells found. Fitting profile by interpolation with the surrounding sane cells...")
		dr = np.diff(r)
		dr_upd = np.append(dr,1.0)
		r = r[dr_upd>0.]
		v = v[dr_upd>0.]
		rho = rho[dr_upd>0.]
		try:
			spline_v = scipl.CubicSpline(r, v)
		except:
			return ValueError("CSM profile generation failed, with excessive radius inversion in CSM cells. Try smaller --tinj (with --skip-eruption to avoid re-running the eruption calculation from scratch).")
	nsize = len(r)
	(rmin, rmax) = (r[0], r[-1])
	# initial guess for outer CSM power-law index (Tsuna & Takei 23, PASJ 75, L19)
	nout_init = 9.

	popt, pcov = curve_fit(CSMprof_func_stripped, np.log(r), np.log(rho), p0 = [1e15, 1e-15, nout_init])
	(r_star, rho_star, nout) = (popt[0], popt[1], popt[2])
	r_remesh = np.logspace(np.log10(rmin*1.001), np.log10(rout*1.001), Ncell)
	rho_remesh = np.zeros(Ncell)
	v_remesh = np.zeros(Ncell)

	for i, x in enumerate(r_remesh):
		if x <= rmax:
			rho_remesh[i] = np.exp(CSMprof_func_stripped(np.log(x), r_star, rho_star, nout))
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
		s = '# abundances\n'
		f.write(s)
		s = 'H  {:.4e}\nHe {:.4e}\nC  {:.4e}\nO  {:.4e}\n'.format(h1,he4,c12,o16)
		f.write(s)

# return characteristic values from the self-similar solution by Chevalier (1982).
# data[:,2] = A, and data[:,3], data[:,4] are values related to the internal energy stored in the shocked ejecta and CSM.
def interpolate_self_similar_solution(n, s):
	data = np.loadtxt("./input/LC/integral.txt")
	x, y, A, Sr, Sf = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4]

	if n < min(x):
		n = min(x)
	elif n > max(x):
		n = max(x)
	
	if s < min(y):
		s = min(y)
	elif s > max(y):
		s = max(y)

	A1  = griddata((x, y), A,  (n, s), method='cubic')
	Sr1 = griddata((x, y), Sr, (n, s), method='cubic')
	Sf1 = griddata((x, y), Sf, (n, s), method='cubic')

	result = np.array([[A1, Sr1, Sf1]])
	np.savetxt('LCFiles/interp_self_similar_values.txt', result, fmt="%1.4e")
