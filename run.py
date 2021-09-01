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
from input.TOPS import gen_op_tbl 
from input.TOPS_multigroup import gen_op_frq
import lightcurve


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
	parser.add_option("--Eexp", metavar = "float", type = "float", action = "append", help = "Explosion energy in erg. This option can be given multiple times (default: 1e51, 3e51, 1e52).")
	parser.add_option("--eruption-innerMr", metavar = "float", type = "float", default=-1.0, help = "The innermost mass coordinate where the energy is injected. If a negative value is given, it sets by default to just outside the helium core.")
	parser.add_option("--inlist-file", metavar = "filename", help = "Inlist file with the ZAMS mass and metallicity information.")
	parser.add_option("--skip-mesa", action = "store_true", help = "Use stellar models pre-computed for input to the mass eruption code.")
	parser.add_option("--analytical-CSM", action = "store_true", default=False, help = "Calibrate CSM by analytical profile given in Tsuna et al (2021). The adiabatic CSM profile is extrapolated to the inner region, correcting the profile obtained from adiabatic calculation that includes artificial shock-compression.")
	parser.add_option("--steady-wind", metavar = "string", type = "string", default='RSGwind', help = "Set how the steady wind CSM is attached to the erupted material. Must be 'attach' or 'RSGwind'.")
	parser.add_option("--calc-multiband", action = "store_true", default=False, help = "Additionally conduct ray-tracing calculations to obtain multi-band light curves (default: false). This calculation is computationally much heavier than obtaining just the bolometric light curve.")

	options, filenames = parser.parse_args()
	available_masses = [13.,14.,15.,16.,17.,18.,19.,20.,22.,24.,26.]
	available_mesa_models = [(mass, 1.) for mass in available_masses]
	if options.zams_m <= 0 or options.zams_z <= 0 or options.tinj <= 0 or options.finj <= 0:
		raise ValueError("The input parameters must all be positive.")
	if options.skip_mesa and (options.zams_m,options.zams_z) not in available_mesa_models:
		print("(M, Z) = (%.1f Msun, %.1f Zsun) not available. Running mesa calculation instead..." % (options.zams_m,options.zams_z))
		options.skip_mesa = False
	# set default value if explosion energy is empty
	if not options.Eexp:
		options.Eexp = [1e51, 3e51, 1e52]

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
convert.convertForHydro(file_me, file_hydro, options.eruption_innerMr)


# energy injection timescale
injectDuration = 1e3 # unit in second
# continueTransfer can be set to true, if radiative transfer scheme needs to be continued even after the eruption.
# However, the computation will be much slower.
convert.setSnhydParam(hydroNumMesh, time_CSM, inject_duration, options.finj, continueTransfer=False)


# run eruptive mass-loss rad-hydro calculation
subprocess.call("./eruption", stdout=open(os.devnull,'wb'))


# obtain light curve at mass eruption
mass_eruption_lc_file = 'LCFiles/mass_eruption_lightcurve.txt'
utils.get_mass_eruption_lightcurve(mass_eruption_lc_file)


#################################################################
#								#
#		IIn light curve model of TS20			#
#								#
#################################################################

# outer extent of the CSN to feed into the LC calculation
r_out = 3e16
CSM_file = 'LCFiles/CSM.txt'
profile_at_cc = 'EruptionFiles/atCCSN.txt'
Y_He, r_edge = utils.remesh_CSM(r_out, profile_at_cc, CSM_file, file_me, analytical_CSM = options.analytical_CSM, steady_wind=options.steady_wind)

# extract the ejecta parameters
Mej, n, delta = utils.calculate_ejecta(file_cc, profile_at_cc, CSM_file, r_edge)

# obtain opacity 
opacity_file = 'LCFiles/opacity.txt'
gen_op_tbl.gen_op_tbl_sct(Y_He, opacity_file)
opacity_file = 'LCFiles/kappa_p.txt'
gen_op_tbl.gen_op_tbl_abs(Y_He, opacity_file)
op_freq_dir = 'LCFiles/opacity_frq'
subprocess.call(["rm", "-rf", op_freq_dir])
subprocess.call(["mkdir", op_freq_dir])
gen_op_frq.gen_op_frq(Y_He, op_freq_dir)


for Eexp in options.Eexp:
	# luminosity at shock
	dir_name_shockprofiles = "LCFiles/ShockProfilesandSpecFiles_"+str(Eexp)
	subprocess.call(["rm", "-r", dir_name_shockprofiles])
	subprocess.call(["mkdir", dir_name_shockprofiles])
	shock_file = 'LCFiles/shock_output_'+str(Eexp)+'erg.txt'
	lightcurve.shock(Eexp, Mej*1.99e33, n, delta, CSM_file, shock_file, dir_name_shockprofiles)

	# radiation transfer
	# bolometric light curve
	IIn_lc_file = 'LCFiles/IIn_lightcurve_'+str(Eexp)+'erg.txt'
	# multi-band light curve if requested
	if options.calc_multiband:
		IIn_lc_band_file = 'LCFiles/IIn_lightcurve_'+str(Eexp)+'erg_mag.txt'
	else:
		IIn_lc_band_file = ''
	lightcurve.transfer(Eexp, Mej*1.99e33, n, delta, r_out, CSM_file, shock_file, IIn_lc_file, IIn_lc_band_file, dir_name_shockprofiles)

	# obtain peak luminosity and rise/decay time in days
	# the rise (decay) time is defined by between peak time and the time when the luminosity first rises(decays) to 1% of the peak.
	try:
		utils.extract_peak_and_rise_time(IIn_lc_file, frac=0.01)
	except ValueError:
		pass
