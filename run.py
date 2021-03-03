from __future__ import print_function
import fileinput
import glob
from optparse import OptionParser
import os
import subprocess
import sys

# our modules
import convert
import utils
import lightcurve
from TOPS import gen_op_tbl 


def parse_command_line():
	parser = OptionParser(
		description = '''Execution script. The MESA calculation will be conducted by running the code with ZAMS mass/Z and inlist file as arguments, For example, to calculate the evolution of a 15Msun star with metallicity 0.02 (~ solar), with inlist file containing these parameters, run the following:\n
		python run.py --zams-m 15 --zams-z 0.02 --inlist-file /path/to/inlist_file
		'''
	)
	parser.add_option("--zams-m", metavar = "float", type = "float", help = "Initial mass in units of solar mass.")
	parser.add_option("--zams-z", metavar = "float", type = "float", help = "Initial metallicity in units of solar metallicity (Z=0.014).")
	parser.add_option("--delta-t", metavar = "float", type = "float", help = "Time from mass eruption to core-collapse, in units of years.")
	parser.add_option("--inlist-file", metavar = "filename", help = "Inlist file with the ZAMS mass and metallicity information.")
	parser.add_option("--skip-mesa", action = "store_true", help = "Use stellar models pre-computed for input to the mass eruption code.")

	options, filenames = parser.parse_args()
	available_masses = [13.,14.,15.,16.,17.,18.,19.,20.,22.,24.,26.,28.,30.]
	available_mesa_models = [(mass, 1.) for mass in available_masses]
	if options.skip_mesa and (options.zams_m,options.zams_z) not in available_mesa_models:
		print("(M, Z) = (%.1f Msun, %.1f Zsun) not available. Running mesa calculation instead..." % (options.zams_m,options.zams_z))
		options.skip_mesa = False

	return options, filenames

# get command line arguments
options, filenames = parse_command_line()


if options.skip_mesa:
	file_cc = 'mesa_models/'+str(int(options.zams_m))+'Msun_Z'+str(0.014*options.zams_z)+'_preccsn.data'
	file_me = file_cc 
	time_CSM = options.delta_t
else:
	#################################################################
	#								#
	#			MESA calculation			#
	#								#
	#################################################################

	# edit the file with the given input zams mass and metallicity
	# FIXME this doesn't work for inlist file without these parameters set.
	for line in fileinput.input(options.inlist_file, inplace=1):
		if 'initial_mass' in line:
			print("      initial_mass = %f" % float(options.zams_m))
		elif 'initial_z' in line:
			print("      initial_z = %f" % (0.014*float(options.zams_z)))
		elif 'Zbase' in line:
			print("      Zbase = %f" % (0.014*float(options.zams_z)))
		else:
			print(line.rstrip())

	# compile mesa script
	mesa_dir = os.path.dirname(options.inlist_file)
	os.chdir(mesa_dir)
	subprocess.call("./mk")

	# run mesa script
	subprocess.call("./rn")
	os.chdir("../")

	# find data file at mass eruption and core collapse. 
	# FIXME we set the mass eruption to 5 years before collapse
	file_cc = mesa_dir+'/pre_ccsn.data'
	file_me = utils.find_mass_eruption(glob.glob(mesa_dir+'/LOGS_to_si_burn/profile*.data'), file_cc, options.delta_t)

	# obtain the time from end of rad-hydro calculation to core-collapse (in years)
	time_CSM = utils.get_mass_eruption_to_core_collapse(file_me, file_cc)
print("from mass eruption to core collapse: %e yrs" % time_CSM, file=sys.stderr)

#################################################################
#                                                               #
#               Eruptive mass loss model of KS20                #
#                                                               #
#################################################################



# convert data for hydro in KS20
file_hydro = 'InitForHydro.txt'
hydroNumMesh = 10000
logscaleRemesh = False

massCutByHand = False # If true, massCutPoint is used. If false, helium core is cutted automatically.
massCutPoint = 1.3 # unit in Msun

subprocess.call(["rm", "f/inclmn.f"])
subprocess.call(["rm", "f/eruptPara.d"])
subprocess.call(["rm", "InitForHydro.txt"])

convert.convertForHydro(file_me, file_hydro, hydroNumMesh, massCutByHand, massCutPoint, logscaleRemesh)



# decide energy injection timescale and luminosity
# find when to calculate the eruptive mass-loss
injectedEnergy = 1.5e47 # unit in erg
injectDuration = 1e3 # unit in second

ScaledByEnvelopeEnergy = True # If enabled, (-1)*(injectedEnergyRate)*(total energy of the envelope) is deposited instead of injectedEnergy.
injectedEnergyRate = 0.3 # around 0.3 is recommended

convert.setSnhydParam(hydroNumMesh,time_CSM,injectedEnergy,injectDuration, ScaledByEnvelopeEnergy, injectedEnergyRate)


# compile eruptive mass-loss rad-hydro calculation (It will be modified to use gfortran later. (Comment by Kuriyama))
subprocess.call(["mkdir", "-p", "snhydOutput"])
subprocess.call(["make", "clean"])
subprocess.call("make")



# run eruptive mass-loss rad-hydro calculation
subprocess.call("./runsnhyd")


#################################################################
#								#
#		IIn light curve model of TS20			#
#								#
#################################################################

# outer extent of the CSN to feed into the LC calculation
r_out = 9.9e15
CSM_file = 'inp-data/CSM.txt'
Y_He, r_edge = utils.remesh_CSM(r_out, 'snhydOutput/atCCSN.txt', CSM_file, file_me)

# extract the ejecta parameters
Mej, n, delta = utils.calculate_ejecta(file_cc, 'snhydOutput/result99.txt', r_edge)
Eexp = 1e51

# obtain opacity 
opacity_file = 'inp-data/opacity.txt'
gen_op_tbl.gen_op_tbl(Y_He, opacity_file)

# luminosity at shock
shock_file = 'inp-data/shock_output.txt'
lightcurve.shock(Eexp, Mej*1.99e33, n, delta, CSM_file, shock_file)

# radiation transfer
IIn_lc_file = 'outp-data/IIn_lightcurve.txt'
lightcurve.transfer(r_out, CSM_file, shock_file, IIn_lc_file)

# obtain peak luminosity and rise/decay time in days
# the rise (decay) time is defined by between peak time and the time when the luminosity first rises(decays) to 1% of the peak.
utils.extract_peak_and_rise_time(IIn_lc_file, frac=0.01)

# obtain light curve at mass eruption
mass_eruption_lc_file = 'outp-data/mass_eruption_lightcurve.txt'
utils.get_mass_eruption_lightcurve(mass_eruption_lc_file)
