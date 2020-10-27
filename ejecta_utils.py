# utilities for extracting ejecta parameters

import math
import mesa_reader as mr
import numpy as np
import warnings


MSUN = 1.9884e+33
RSUN = 6.96e+10
G = 6.6743e-8

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
		lgpenv = [lgp for i, lgp in enumerate(data.logp) if data.mass[i] > He_core_mass]
	# otherwise, it should be the helium envelope
	else:
		lgrhoenv = [lgrho for i, lgrho in enumerate(data.logRho) if data.mass[i] > CO_core_mass]
		lgpenv = [lgp for i, lgp in enumerate(data.logp) if data.mass[i] > CO_core_mass]
	return total_mass, CO_core_mass, lgrho_env, lgp_env


# ejecta calculation script
def calculate_ej_from_mesa(data_file):
	# remnant mass
	total_mass, CO_core_mass, lgrho_env, lgp_env = cc_param_extractor(data_file)
	remnant_mass = remnant_from_CO(CO_core_mass)
	Mej = total_mass - remnant_mass
	# obtain n from fitting rho vs p at envelope. We fit with P = K*rho^N, where N is the polytripic index.
	Npol = np.polyfit(lgrho_env, lgp_env)[0]
	# Determine n using Matzner & McKee 1999 eq 25
	beta = 0.19
	n = (Npol+1.+3.*beta*Npol)/(beta*Npol)
	# set delta to 1 for now
	delta = 1.0
	return Mej, n, delta


# obtain time from mass eruption till core-collapse
def get_end_of_sim_to_core_collapse(data_file_at_mass_eruption, data_file_at_core_collapse):
	data_me = mr.MesaData(data_file_at_mass_eruption)
	data_cc = mr.MesaData(data_file_at_core_collapse)
	# in years
	return data_cc.star_age - data_me.star_age


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
	rmin = r_in[0]
	X_edge = CSM[-1,5]
	Y_edge = CSM[-1,6]
	vwind = 1.6 *  math.sqrt(2.*G*data.star_mass*MSUN/data.photosphere_r/RSUN)
	wind_Mdot_vw = -data.star_mdot * MSUN / 3.15e7 / vwind 

	rs = np.logspace(math.log10(rmin*1.001), math.log10(rmax*1.001), 1000)
	for i, r in enumerate(rs):
		if r < r_in[-1]:
			# obtain rho by log interpolation, Mr and v by linear. Mr is not used anyway
			index = len([thisr for thisr in r_in if thisr < r])-1
			fraction = (r - r_in[index]) / (r_in[index+1] - r_in[index])
			rho = CSM[index, 4]**(1.-fraction) * CSM[index+1,4]**(fraction)
			Mr = CSM[index, 1]*(1.-fraction) + CSM[index+1,1]*(fraction)
			v = CSM[index, 3]*(1.-fraction) + CSM[index+1,3]*(fraction)
			X = CSM[index, 5]*(1.-fraction) + CSM[index+1,5]*(fraction)
			Y = CSM[index, 6]*(1.-fraction) + CSM[index+1,6]*(fraction)
			print("%d %.4g %.4g %.4g %.4g %.4g %.4g" % (i, Mr, r, v, rho, X, Y), file=fout)
		else:
			# use the wind profile at that point
			rho = wind_Mdot_vw / (4.*math.pi*r**2)
			Mr += 4.*math.pi*r**2*rho*(r-rs[i-1])
			# X, Y are the edge value
			print("%d %.4g %.4g %.4g %.4g %.4g %.4g" % (i, Mr, r, vwind, rho, X_edge, Y_edge), file=fout)
	fout.close()

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
	return peakL, risetime, decaytime 
