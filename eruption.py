from __future__ import print_function
import fileinput
import glob
from optparse import OptionParser
import os
import subprocess
import sys

# our modules
from utils import convert
from utils import utils

def parse_command_line():
	parser = OptionParser(
		description = '''Execution script. The MESA calculation will be conducted by running the code with ZAMS mass/Z and inlist file as arguments, For example, to calculate the evolution of a 15Msun star with metallicity 0.02 (~ solar), with inlist file containing these parameters, run the following:\n
		python run.py --zams-m 15 --zams-z 0.02 --inlist-file /path/to/inlist_file
		'''
	)
	parser.add_option("--zams-m", metavar = "float", type = "float", help = "Initial mass in units of solar mass.")
	parser.add_option("--zams-z", metavar = "float", type = "float", help = "Initial metallicity in units of solar metallicity (Z=0.014).")
	parser.add_option("--tinj", metavar = "float", type = "float", help = "Time from mass eruption to core-collapse, in units of years.")
	parser.add_option("--finj", metavar = "float", type = "float", default=0.3, help = "Energy injected at the base of the stellar envelope, scaled with the envelope's binding energy (default: 0.3).")
	parser.add_option("--inject-duration", metavar = "float", type = "float", default=1e3, help = "Duration of the energy injection at the base of the stellar envelope (default: 1000 seconds, a value much smaller than the envelope's dynamical timescale).")
	parser.add_option("--inlist-file", metavar = "filename", help = "Inlist file with the ZAMS mass and metallicity information.")
	parser.add_option("--skip-mesa", action = "store_true", help = "Use stellar models pre-computed for input to the mass eruption code.")

	options, filenames = parser.parse_args()
	available_masses = [13.,14.,15.,16.,17.,18.,19.,20.,22.,24.,26.]
	available_mesa_models = [(mass, 1.) for mass in available_masses]
	if options.zams_m <= 0 or options.zams_z <= 0 or options.tinj <= 0 or options.finj <= 0:
		raise ValueError("The input parameters must all be positive.")
	if options.skip_mesa and (options.zams_m,options.zams_z) not in available_mesa_models:
		print("(M, Z) = (%.1f Msun, %.1f Zsun) not available. Running mesa calculation instead..." % (options.zams_m,options.zams_z))
		options.skip_mesa = False

	return options, filenames

# get command line arguments
options, filenames = parse_command_line()


if options.skip_mesa:
	file_cc = 'input/mesa_models/'+str(int(options.zams_m))+'Msun_Z'+str(0.014*options.zams_z)+'_preccsn.data'
	file_me = file_cc 
	time_CSM = options.tinj
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
	file_me = utils.find_mass_eruption(glob.glob(mesa_dir+'/LOGS_to_si_burn/profile*.data'), file_cc, options.tinj)

	# obtain the time from end of rad-hydro calculation to core-collapse (in years)
	time_CSM = utils.get_mass_eruption_to_core_collapse(file_me, file_cc)
print("from mass eruption to core collapse: %e yrs" % time_CSM, file=sys.stderr)

#################################################################
#                                                               #
#               Eruptive mass loss model of KS20                #
#                                                               #
#################################################################



# convert data for hydro in KS20
file_hydro = 'EruptionFiles/InitForHydro.txt'
hydroNumMesh = 10000
logscaleRemesh = False

massCutByHand = False # If true, massCutPoint is used. If false, helium core is cutted automatically.
massCutPoint = 1.3 # unit in Msun

subprocess.call(["rm", "src/eruption/f/inclmn.f"])
subprocess.call(["rm", "src/eruption/f/eruptPara.d"])
subprocess.call(["rm", "EruptionFiles/InitForHydro.txt"])

convert.convertForHydro(file_me, file_hydro, hydroNumMesh, massCutByHand, massCutPoint, logscaleRemesh)


# decide energy injection timescale and luminosity
injectDuration = options.inject_duration # unit in second
injectedEnergyRate = options.finj # around 0.3 is recommended

# Stop radiative transfer calculation from well after mass eruption.
# If true, radiative transfer scheme is activated even after the eruption.
continueTransfer = False
# FIXME remove extra argument of 0, True
convert.setSnhydParam(hydroNumMesh, time_CSM, 0, injectDuration, True, injectedEnergyRate, continueTransfer)


# compile eruptive mass-loss rad-hydro calculation (It will be modified to use gfortran later. (Comment by Kuriyama))
os.chdir("src/eruption")
subprocess.call(["make", "clean"])
subprocess.call("make")
os.chdir("../../")


# run eruptive mass-loss rad-hydro calculation
subprocess.call("./eruption", stdout=open(os.devnull,'wb'))


# obtain light curve at mass eruption
mass_eruption_lc_file = 'LCFiles/mass_eruption_lightcurve.txt'
utils.get_mass_eruption_lightcurve(mass_eruption_lc_file)
