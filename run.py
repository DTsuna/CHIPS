from __future__ import print_function
import fileinput
from optparse import OptionParser
import subprocess

def parse_command_line():
	parser = OptionParser(
		description = '''Execution script. The MESA calculation will be conducted by running the code with ZAMS mass/Z and inlist file as arguments, For example, to calculate the evolution of a 15Msun star with metallicity 0.02 (~ solar), with inlist file containing these parameters in this_inlist_file, run the following:\n
		python run.py --zams-m 15 --zams-z 0.02 --inlist-file this_inlist_file
		'''
	)
	parser.add_option("--zams-m", metavar = "float", type = "float", help = "Initial mass.")
	parser.add_option("--zams-z", metavar = "float", type = "float", help = "Initial metallicity.")
	parser.add_option("--inlist-file", metavar = "filename", help = "Inlist file with the ZAMS mass and metallicity information.")
	options, filenames = parser.parse_args()

	return options, filenames

# get command line arguments
options, filenames = parse_command_line()


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
		print("      initial_z = %f" % float(options.zams_z))
	elif 'Zbase' in line:
		print("      Zbase = %f" % float(options.zams_z))
	else:
		print(line.rstrip())

# compile mesa script
#subprocess.call("./mk")

# run mesa script
#subprocess.call("./rn")



#################################################################
#								#
#		Eruptive mass loss model of KS20		#
#								#
#################################################################

# find when to calculate the eruptive mass-loss

# decide energy injection timescale and luminosity

# do eruptive mass-loss rad-hydro calculation


#################################################################
#								#
#		Calculate CSM at core-collapse			#
#								#
#################################################################

# obtain the time from end of rad-hydro calculation to core-collapse

# read in CSM data

# orbit calculation w/ gravity of central star

# record CSM data at core-collapse



#################################################################
#								#
#		IIn light curve model of TS20			#
#								#
#################################################################

# obtain the ejecta & CSM parameters

# do light curve calculation

# record light curve

# obtain peak luminosity and timescales
