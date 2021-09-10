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
	parser.add_option("--eruption-innerMr", metavar = "float", type = "float", default=-1.0, help = "The innermost mass coordinate where the energy is injected. If a negative value is given, it sets by default to just outside the helium core.")
	parser.add_option("--inject-duration", metavar = "float", type = "float", default=1e3, help = "Duration of the energy injection at the base of the stellar envelope (default: 1000 seconds, a value much smaller than the envelope's dynamical timescale).")
	parser.add_option("--inlist-file", metavar = "filename", help = "Inlist file with the ZAMS mass and metallicity information.")
	parser.add_option("--skip-mesa", action = "store_true", help = "Use stellar models pre-computed for input to the mass eruption code.")

	options, filenames = parser.parse_args()
	available_masses = [13.,14.,15.,16.,17.,18.,19.,20.,22.,24.,26.]
	available_mesa_models = [(mass, 1.) for mass in available_masses]
	if options.zams_m <= 0 or options.zams_z <= 0 or options.tinj <= 0 or options.finj <= 0:
		raise ValueError("The input parameters zams_m, zams_z, tinj, finj must all be positive.")
	if options.skip_mesa and (options.zams_m,options.zams_z) not in available_mesa_models:
		print("(M, Z) = (%.1f Msun, %.1f Zsun) not available. Running mesa calculation instead..." % (options.zams_m,options.zams_z))
		options.skip_mesa = False

	return options, filenames

# get command line arguments
options, filenames = parse_command_line()

if options.skip_mesa:
	file_cc = 'input/mesa_models/'+str(int(options.zams_m))+'Msun_Z'+str(0.014*options.zams_z)+'_preccsn.data'
else:
	#################################################################
	#								#
	#			MESA calculation			#
	#								#
	#################################################################

	# edit the file with the given input zams mass and metallicity
	# NOTE the inlist file needs to already contain "initial_mass", "initial_z" and "Zbase", and "filename_for_profile_when_terminate".
	for line in fileinput.input(options.inlist_file, inplace=1):
		initial_mass_flag = 0
		initial_z_flag = 0
		Zbase_flag = 0
		ccfile_flag = 0
		if 'initial_mass' in line:
			initial_mass_flag = 1
			print("      initial_mass = %f" % float(options.zams_m))
		elif 'initial_z' in line:
			initial_z_flag = 1
			print("      initial_z = %f" % (0.014*float(options.zams_z)))
		elif 'Zbase' in line:
			Zbase_flag = 1
			print("      Zbase = %f" % (0.014*float(options.zams_z)))
		elif 'filename_for_profile_when_terminate' in line:
			ccfile_flag = 1
			# extract filename for pre-ccSN model
			preSN_filename = line.split()[-1].split('=')[-1].strip('\'') 
			print(line.rstrip())
		else:
			print(line.rstrip())
	
	assert (initial_mass_flag, initial_z_flag, Zbase_flag, ccfile_flag) == (1, 1, 1, 1), "arguments (initial_mass, initial_z, Zbase, filename_for_profile_when_terminate) all have to be included in the inlist file" 

	# compile mesa script
	mesa_dir = os.path.dirname(options.inlist_file)
	os.chdir(mesa_dir)
	subprocess.call("./mk")

	# run mesa script
	subprocess.call("./rn")
	os.chdir("../")

	# find data file at mass eruption, assumed to be same as core collapse. 
	file_cc = mesa_dir + '/' + preSN_filename


#################################################################
#                                                               #
#               Eruptive mass loss model of KS20                #
#                                                               #
#################################################################


# convert data for hydro in KS20
file_eruption = 'EruptionFiles/InitForHydro.txt'
convert.convertForEruption(file_cc, file_eruption, massCutPoint=options.eruption_innerMr)

# continueTransfer can be set to true, if radiative transfer scheme needs to be continued even after the eruption.
# However, the computation will be much slower.
convert.setEruptionParam(options.tinj, options.inject_duration, options.finj, continueTransfer=False)

# run eruptive mass-loss rad-hydro calculation
subprocess.call("./eruption", stdout=open(os.devnull,'wb'))

# obtain light curve at mass eruption
mass_eruption_lc_file = 'LCFiles/mass_eruption_lightcurve.txt'
utils.get_mass_eruption_lightcurve(mass_eruption_lc_file)
