from __future__ import print_function
import fileinput
from optparse import OptionParser
import subprocess

import ejecta_utils
import convert


def parse_command_line():
	parser = OptionParser(
		description = '''Execution script. The MESA calculation will be conducted by running the code with ZAMS mass/Z and inlist file as arguments, For example, to calculate the evolution of a 15Msun star with metallicity 0.02 (~ solar), with inlist file containing these parameters, run the following:\n
		python run.py --zams-m 15 --zams-z 0.02 --inlist-file /path/to/inlist_file
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

# find data file at core collapse
file_cc = 'pre_ccsn.data'
#file_cc = 'profile54.data'

# obtain the time from end of rad-hydro calculation to core-collapse (in years)
time_CSM = ejecta_utils.get_mass_eruption_to_core_collapse(file_me, file_cc)

#################################################################
#                                                               #
#               Eruptive mass loss model of KS20                #
#                                                               #
#################################################################



# convert data for hydro in KS20
file_hydro = 'InitForHydro.txt'
hydroNumMesh = 10000

subprocess.call(["rm", "snhyd/inclmn.f"])
subprocess.call(["rm", "snhyd/eruptPara.d"])
subprocess.call(["rm", "InitForHydro.txt"])

convert.convertForHydro(file_cc, file_hydro, hydroNumMesh)



# decide energy injection timescale and luminosity
# find when to calculate the eruptive mass-loss
injectedEnergy = 1.5e47
injectDuration = 1e3

convert.setSnhydParam(hydroNumMesh,time_CSM,injectedEnergy,injectDuration)



# compile eruptive mass-loss rad-hydro calculation (It will be modified to use gfortran later. (Comment by Kuriyama))
subprocess.call(["mkdir", "snhydOutput"])
subprocess.call(["make", "clean"])
compileEOS = "gfortran -fno-automatic  -c -o snhyd/eos_helm.o snhyd/eos_helm.f90"
subprocess.call(compileEOS.split())
compileEOS = "gfortran -fno-automatic  -c -o snhyd/eos_helm_e.o snhyd/eos_helm_e.f90"
subprocess.call(compileEOS.split())
compileEOS = "gfortran -fno-automatic  -c -o snhyd/eos_helm_p.o snhyd/eos_helm_p.f90"
subprocess.call(compileEOS.split())
subprocess.call("make")



# run eruptive mass-loss rad-hydro calculation
subprocess.call("./runsnhyd")



#################################################################
#								#
#		Calculate CSM at core-collapse			#
#								#
#################################################################

# CSM distribution will be calculated by hydro until core-collapse

#################################################################
#								#
#		IIn light curve model of TS20			#
#								#
#################################################################

# extract the ejecta parameters
Mej, n, delta = ejecta_utils.calculate_ej_from_mesa(file_cc)
explosionEnergy = 1e51

# do light curve calculation

# record light curve

# obtain peak luminosity and rise/decay time in days
# the rise (decay) time is defined by between peak time and the time when the luminosity first rises(decays) to 1% of the peak.
peakL, riset, decayt = ejecta_utils.extract_peak_and_rise_time(LC_output, frac=0.01)
