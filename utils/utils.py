# utilities functions for run.py 

import math
import mesa_reader as mr
import numpy as np
import sys
import warnings
import os
import re
# check of whether scipy can be imported
try:
	from scipy.optimize import curve_fit
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


MSUN = 1.9884e+33
RSUN = 6.96e+10
G = 6.6743e-8

# for CSM density profile fitting
nmax = 12.
def CSMprof_func(r, r_break, rho_break, yrho):
	return np.log( rho_break * (( (np.exp(r) / r_break)**(1.5/yrho) + (np.exp(r) / r_break)**(nmax/yrho) ) /2. )**(-yrho) )

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
		Mrem = 10.**(-0.02466*CO_core_mass+1.28070)
	return Mrem


# extractor of envelope profile from the script 
def cc_param_extractor(data_file):
	data = mr.MesaData(data_file)
	total_mass = data.star_mass
	CO_core_mass = data.c_core_mass
	return total_mass, CO_core_mass


# ejecta calculation script
# r_edge is the edge of the ejecta, i.e. start of the CSM
def calculate_ejecta(data_file, file_at_cc, file_CSM, r_edge):
	renv = np.loadtxt(file_at_cc, skiprows=1)[:,2]
	rhoenv = np.loadtxt(file_at_cc, skiprows=1)[:,4]
	penv =  np.loadtxt(file_at_cc, skiprows=1)[:,7]
	lgrhoenv = [math.log(rho) for i, rho in enumerate(rhoenv) if renv[i]<r_edge]
	lgpenv = [math.log(p) for i, p in enumerate(penv) if renv[i]<r_edge]
	# obtain polytropic index from fitting rho vs p at envelope. We fit with P = K*rho^(1+1/N), where N is the polytripic index.
	gamma = np.polyfit(lgrhoenv, lgpenv, 1)[0]
	Npol = 1./(gamma-1.)
	# Determine n using Matzner & McKee 1999 eq 25
	beta = 0.19
	n = (Npol+1.+3.*beta*Npol)/(beta*Npol)
	assert n>5, "Ejecta index should be >5, now %e" % n
	# set delta to 1 for now
	delta = 1.0
	# ejecta mass
	total_mass, CO_core_mass = cc_param_extractor(data_file)
	remnant_mass = remnant_from_CO(CO_core_mass)
	# obtain CSM mass from CSM file
	MrCSM = np.loadtxt(file_CSM, skiprows=1)[:,1]
	CSM_mass = (MrCSM[-1] - MrCSM[0]) / MSUN 
	Mej = total_mass - remnant_mass - CSM_mass
	assert Mej > 0.0
	print("Mej:%f Msun, n:%f, delta:%f" % (Mej, n, delta), file=sys.stderr)
	return Mej, n, delta


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
	outward_global_slope = [r_in[i+10]/rho_in[i+20]*(rho_in[i+20]-rho_in[i+10])/(r_in[i+20]-r_in[i+10]) for i,r in enumerate(r_in) if i<len(r_in)-20]
	istop = max([i for i in range(len(slope[:-21])) if slope[i]<-100 and outward_global_slope[i+1]>-4.0 and outward_global_slope[i+1]<-1.0])
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
			# we fix the density profile as rho\propto r^(-1.5) inside the radius where the CSM density become
			# unreliable due to artificial shocks. this profile is merely a guess: but should be somewhat accurate
			# well close to the stellar surface
			if analytical_CSM and scipy_exists:
				rho_out[i] = math.exp(CSMprof_func(math.log(r), r_break, rho_break, yrho))
			else:
				rho_out[i] = rho_in[istop] * (r/rstop)**(-1.5)
			# Mr is obtained later once rho and r are determined
			# for other ones use values at istop, since they are not important anyway
			v_out[i] = CSM[istop, 3]
			X_out[i] = CSM[istop, 5]
			Y_out[i] = CSM[istop, 6]
			counter += 1
		elif r < r_in[-1]:
			# initial Mr is that at istop.
			# NOTE this sets the zero point of the enclosed mass, which is presumably close
			# to the actual value. We don't really mind because we only use the difference of Mr
			# between 2 cells, but if we include the central star's gravity for calculation of
			# the light curve, this has to be revisited.
			last_Mr = CSM[istop, 1]
			# obtain rho by log interpolation, Mr and v by linear. Mr is not used anyway
			index = len([thisr for thisr in r_in if thisr < r])-1
			fraction = (r - r_in[index]) / (r_in[index+1] - r_in[index])
			v_out[i] = CSM[index, 3]*(1.-fraction) + CSM[index+1,3]*(fraction)
			rho_out[i] = CSM[index, 4]**(1.-fraction) * CSM[index+1,4]**(fraction) 
			Mr_out[i] = last_Mr + 4.*math.pi*r**2*rho_out[i]*(r-r_out[i-1])
			X_out[i] = CSM[index, 5]*(1.-fraction) + CSM[index+1,5]*(fraction)
			Y_out[i] = CSM[index, 6]*(1.-fraction) + CSM[index+1,6]*(fraction)
			# for next step
			Y_avrg = (Y_avrg*last_Mr + Y_out[i]*(Mr_out[i]-last_Mr))/Mr_out[i]
			last_Mr = Mr_out[i]
		else:
			# use a profile that connects to the wind profile with a Gaussian drop 
			if steady_wind == 'RSGwind':
				rho_out[i] = rho_in[-1]*math.exp(1.-(r/r_in[-1])**2)  + wind_Mdot_vw / (4.*math.pi) * (1./r**2)
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
			Y_avrg = (Y_avrg*last_Mr + Y_out[i]*(Mr_out[i]-last_Mr))/Mr_out[i]
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

	# extract Y_avrg, needed for opacity calculation, and start of CSM, needed to set end of ejecta for ejecta calculation
	return Y_avrg, max(r_out[0], rstop)

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
	decaytime = time[np.argmin([abs(L-peakL*frac) for i,L in enumerate(lum) if time[i]>peaktime])] - peaktime
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
