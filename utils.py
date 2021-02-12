# utilities functions for run.py 

import math
import mesa_reader as mr
import numpy as np
import sys
import warnings


MSUN = 1.9884e+33
RSUN = 6.96e+10
G = 6.6743e-8


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
	else:
		# FIXME currently extrapolating the NS value to BH mass range, with maximum at 2.1Msun.
		warnings.warn("This CO core mass is predicted to lead to BH formation. Extrapolating the NS relation up to 2.1Msun...")
		Mrem = max(2.1, 0.03357 * CO_core_mass + 1.31780)
	return Mrem


# extractor of envelope profile from the script 
def cc_param_extractor(data_file):
	data = mr.MesaData(data_file)
	total_mass = data.star_mass
	He_core_mass = data.he_core_mass
	CO_core_mass = data.c_core_mass
	# FIXME this criterion should be revisited
	# if hydrogen envelope significantly exists, the env should be the hydrogen envelope
	if He_core_mass < 0.99*total_mass:
		lgrhoenv = [lgrho for i, lgrho in enumerate(data.logRho) if data.mass[i] > He_core_mass]
		lgpenv = [lgp for i, lgp in enumerate(data.logP) if data.mass[i] > He_core_mass]
	# otherwise, it should be the helium envelope
	else:
		lgrhoenv = [lgrho for i, lgrho in enumerate(data.logRho) if data.mass[i] > CO_core_mass]
		lgpenv = [lgp for i, lgp in enumerate(data.logP) if data.mass[i] > CO_core_mass]
	return total_mass, CO_core_mass, lgrhoenv, lgpenv


# ejecta calculation script
def calculate_ej_from_mesa(data_file):
	# remnant mass
	total_mass, CO_core_mass, lgrhoenv, lgpenv = cc_param_extractor(data_file)
	remnant_mass = remnant_from_CO(CO_core_mass)
	Mej = total_mass - remnant_mass
	# obtain polytropic index from fitting rho vs p at envelope. We fit with P = K*rho^N, where N is the polytripic index.
	Npol = np.polyfit(lgrhoenv, lgpenv, 1)[0]
	# Determine n using Matzner & McKee 1999 eq 25
	beta = 0.19
	n = (Npol+1.+3.*beta*Npol)/(beta*Npol)
	# set delta to 1 for now
	delta = 1.0
	print("Mej:%f, n:%f, delta:%f" % (Mej, n, delta), file=sys.stderr)
	return Mej, n, delta


# remesh CSM for input to the light curve calculation
def remesh_CSM(rmax, CSM_in, CSM_out, data_file_at_mass_eruption, Ncell=1000):
	# copy first line
	fout = open(CSM_out, 'w')
	with open(CSM_in, 'r') as fin:
		fout.write(fin.readline())
	
	# record values with edited mesh
	data = mr.MesaData(data_file_at_mass_eruption)
	CSM = np.loadtxt(CSM_in, skiprows=1)
	r_in = CSM[:,2]
	v_in = CSM[:,3]
	rho_in = CSM[:,4]
	rmin = max(r_in[0], 2.*data.photosphere_r*RSUN)
	X_edge = CSM[-1,5]
	Y_edge = CSM[-1,6]
	vwind = 1.6 *  math.sqrt(2.*G*data.star_mass*MSUN/data.photosphere_r/RSUN)
	if data.star_mdot > 0.0:
		wind_Mdot_vw = -data.star_mdot * MSUN / 3.15e7 / vwind 
	else:
		# FIXME input random mass loss rate of 10^(-5)Msun/yr for now
		wind_Mdot_vw = 1e-5 * MSUN / 3.15e7 / vwind 
	
	last_Mr = 0.0
	Y_avrg = 0.0

	rs = np.logspace(math.log10(rmin*1.001), math.log10(rmax*1.001), Ncell)
	# find outermost radius where the velocity suddenly decreases. 
	# this can be the radius where the CSM density becomes unreliable if there exists an artificial shock, or simply
	# can be the boundary of the star and the CSM.
	try:
		# ignore the outermost cell
		istop = max([i for i in range(len(v_in)-1) if v_in[i]-v_in[i-1]<-5e5])
		rstop = r_in[istop]
	except:
		istop = 0
		rstop = 0.0
	for i, r in enumerate(rs):
		if r < rstop:
			# we fix the density profile as rho\propto r^(-1.5) inside the radius where the CSM density become
			# unreliable due to artificial shocks. this profile is merely a guess: but should be somewhat accurate
			# well close to the stellar surface
			rho = rho_in[istop] * (r/rstop)**(-1.5)
			# use values at istop
			# FIXME Mr is incorrect; fix or delete Mr
			Mr = CSM[istop, 1]
			v = CSM[istop, 3]
			X = CSM[istop, 5]
			Y = CSM[istop, 6]
			print("%d %.8g %.8g %.8g %.8g %.8g %.8g" % (i, Mr, r, v, rho, X, Y), file=fout)
		elif r < r_in[-1]:
			# obtain rho by log interpolation, Mr and v by linear. Mr is not used anyway
			index = len([thisr for thisr in r_in if thisr < r])-1
			fraction = (r - r_in[index]) / (r_in[index+1] - r_in[index])
			rho = CSM[index, 4]**(1.-fraction) * CSM[index+1,4]**(fraction)
			Mr = CSM[index, 1]*(1.-fraction) + CSM[index+1,1]*(fraction)
			v = CSM[index, 3]*(1.-fraction) + CSM[index+1,3]*(fraction)
			X = CSM[index, 5]*(1.-fraction) + CSM[index+1,5]*(fraction)
			Y = CSM[index, 6]*(1.-fraction) + CSM[index+1,6]*(fraction)
			Y_avrg = (Y_avrg*last_Mr + Y*(Mr-last_Mr))/Mr
			last_Mr = Mr
			print("%d %.8g %.8g %.8g %.8g %.8g %.8g" % (i, Mr, r, v, rho, X, Y), file=fout)
		else:
			# use a profile that connects to the wind profile with a Gaussian drop 
			#rho = rho_in[-1]*math.exp(1.-(r/r_in[-1])**2)  + wind_Mdot_vw / (4.*math.pi) * (1./r**2)
			# or connect a wind profile to the edge of the erupted CSM profile
			rho = rho_in[-1] * (r_in[-1]/r)**2
			Mr += 4.*math.pi*r**2*rho*(r-rs[i-1])
			Y_avrg = (Y_avrg*last_Mr + Y_edge*(Mr-last_Mr))/Mr
			last_Mr = Mr
			# X, Y are the edge value
			print("%d %.8g %.8g %.8g %.8g %.8g %.8g" % (i, Mr, r, vwind, rho, X_edge, Y_edge), file=fout)
	fout.close()
	return Y_avrg

# extract peak luminosity and rise time, defined as the time from (frac*L_peak) to L_peak
def extract_peak_and_rise_time(LC_file, frac):
	if frac < 0.01:
		return ValueError("Frac too small to give meaningful rise and decay times.")
	time = np.loadtxt(LC_file)[:,0]
	lum = np.loadtxt(LC_file)[:,1]
	peaktime = time[np.argmax(lum)]
	peakL = max(lum)
	risetime = peaktime - time[np.argmin([abs(L-peak*frac) for i,L in enumerate(lum) if time[i]<peaktime])]
	decaytime = time[np.argmin([abs(L-peak*frac) for i,L in enumerate(lum) if time[i]>peaktime])] - peaktime
	print("peak: %e erg/s. rise time: %e days, decay time:%e days" % (peakL, risetime, decaytime), file=sys.stderr)
