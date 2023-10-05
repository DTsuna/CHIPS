import glob
from optparse import OptionParser
import os
import sys
import subprocess

# our modules
from utils import convert
from utils import utils
from input.TOPS import gen_op_tbl 
from input.TOPS_multigroup import gen_op_frq
import lightcurve

MSUN = 1.989E33


def parse_command_line():
	parser = OptionParser(
		description = '''Execution script. The MESA calculation will be conducted by running the code with ZAMS mass/Z and inlist file as arguments, For example, to calculate the interaction-powered supernova of a ZAMS 15Msun, solar metallicity star with tinj=10yrs, finj=0.5, and explosion energy 1e51 ergs, run the following:\n
		python run.py --tinj 10 --finj 0.5 --Eej 1e51 --stellar-model input/mesa_models/15Msun_Z0.014_preccsn.data --analytical-CSM
		'''
	)
	parser.add_option("--tinj", metavar = "float", type = "float", help = "Time from mass eruption to core-collapse, in units of years (required).")
	parser.add_option("--finj", metavar = "float", type = "float", help = "Energy injected at the base of the stellar envelope, scaled with the envelope's binding energy (required).")
	parser.add_option("--Eej", metavar = "float", type = "float", action = "append", help = "Explosion energy in erg. This option can be given multiple times (default: 1e51, 3e51, 1e52).")
	parser.add_option("--injection-duration", metavar = "float", type = "float", default = 1000, help = "Duration of energy injection for mass eruption calculation, in seconds (default: 1000).")
	parser.add_option("--Mni", metavar = "float", type = "float", default = 0., help = "Radioactive nickel 56 mass, in units of solar mass (default: 0 Msun).")
	parser.add_option("--stellar-model", metavar = "filename", help = "Path to the input stellar model (required). This should be one of the stellar model files created after running MESA (which usually end with '.data'.). If --run-mesa is called, this needs to be the stellar model file that you want to provide as input of the CHIPS code (e.g. the file provided by the input 'filename_for_profile_when_terminate' in one of the inlist files.).")
	parser.add_option("--run-mesa", action = "store_true", help = "Call to run MESA in this script and get a new stellar model.")
	parser.add_option("--mesa-path", metavar = "string", type = "string", help = "Path to the execution files of MESA.")
	parser.add_option("--eruption-innerMr", metavar = "float", type = "float", default=-1.0, help = "The innermost mass coordinate where the energy is injected. If no argument or a negative value is given, it sets by default to just outside the helium core.")
	parser.add_option("--analytical-CSM", action = "store_true", default=False, help = "Calibrate CSM by analytical profile given in Tsuna et al (2021). The adiabatic CSM profile is extrapolated to the inner region, correcting the profile obtained from adiabatic calculation that includes artificial shock-compression.")
	parser.add_option("--steady-wind", metavar = "string", type = "string", default='RSGwind', help = "Specify how the steady wind CSM is attached to the erupted material. Must be 'attach' or 'RSGwind' (default: RSGwind). 'attach' simply connects a wind profile to the outermost cell profile, while 'RSGwind' smoothly connects a red supergiant wind to the erupted material.")
	parser.add_option("--calc-multiband", action = "store_true", default=False, help = "Additionally conduct ray-tracing calculations to obtain multi-band light curves (default: false). This calculation is computationally heavier than obtaining just the bolometric light curve.")
	parser.add_option("--opacity-table",  metavar = "filename", help = "A custom opacity table used for the mass eruption calculation. If not called, analytical opacity formula (Kuriyama & Shigeyama 20) is used. The opacity file should have the format like the files in the input/rosseland directory.")

	options, filenames = parser.parse_args()

	# sanity checks
	assert options.stellar_model is not None, "A stellar model file for input to CHIPS has to be provided."
	if options.tinj is None or options.finj is None:
		raise ValueError("Parameters tinj and finj needs to be provided.")
	if options.run_mesa:
		assert options.mesa_path is not None, "A valid existing directory has to be given for --mesa-path if --run-mesa is called."
	if options.tinj <= 0 or options.finj <= 0:
		raise ValueError("The input parameters tinj, finj must both be positive.")
	# set default value if explosion energy is empty
	if not options.Eej:
		options.Eej = [1e51, 3e51, 1e52]

	return options, filenames


# get command line arguments
options, filenames = parse_command_line()

if options.run_mesa:
	#################################################################
	#			MESA calculation			#
	#################################################################

	# compile mesa script
	orig_path = os.getcwd()
	os.chdir(options.mesa_path)
	subprocess.call("./mk")
	# run mesa script
	subprocess.call("./rn")
	os.chdir(orig_path)
	# find desired stellar model file as input to the mass eruption calculation.
	file_cc = options.stellar_model

else:
	file_cc = options.stellar_model


#################################################################
#                                                               #
#               Eruptive mass loss model of KS20                #
#                                                               #
#################################################################

# discriminate progenitor type
# 0: Type IIn
# 1: Type Ibn
# 2: Type Icn
D = utils.discriminantProgModel(file_cc)

# store abundance
utils.genAbundanceTable(file_cc)

# generate tables for mu, T.
lightcurve.opacTable(D)

# convert data for eruption calculation
file_eruption = 'EruptionFiles/InitForHydro.txt'
convert.convertForEruption(file_cc, file_eruption, options.eruption_innerMr, D)

# continueTransfer can be set to true, if radiative transfer scheme needs to be continued even after the eruption. However, the computation will be much slower.
convert.setEruptionParam(options.tinj, options.injection_duration, options.finj, D, continueTransfer=True, OpacityTable=options.opacity_table)

# run eruptive mass-loss rad-hydro calculation
subprocess.call("./eruption")
#subprocess.call("./eruption", stdout=open(os.devnull,'wb'))

# obtain light curve at mass eruption
# comment out the following now that light curve output in EruptionFiles/photosphere.txt
# mass_eruption_lc_file = 'LCFiles/mass_eruption_lightcurve.txt'
# utils.get_mass_eruption_lightcurve(mass_eruption_lc_file)


#################################################################
#								#
#		IIn light curve model of TS20			#
#								#
#################################################################


# outer extent of the CSN to feed into the LC calculation
r_out = 3e16
# remesh CSM in order to correct for shocks in the hydro simulation and extend to r_out.
CSM_file = 'LCFiles/CSM.txt'
profile_at_cc = 'EruptionFiles/atCCSN.txt'
if D == 0:
# obtain opacity
	Y_He = utils.remesh_CSM(r_out, profile_at_cc, CSM_file, file_cc, analytical_CSM = options.analytical_CSM, steady_wind=options.steady_wind)
	opacity_file = 'LCFiles/opacity.txt'
	gen_op_tbl.gen_op_tbl_sct(Y_He, opacity_file)
	opacity_file = 'LCFiles/kappa_p.txt'
	gen_op_tbl.gen_op_tbl_abs(Y_He, opacity_file)
elif D == 1:
	utils.remesh_evolv_CSM(options.tinj, r_out, CSM_file, file_cc, Ncell=1000)
	subprocess.call(["cp", "./input/rosseland/opacity_X0000Y0986Z0014.txt", "./LCFiles/opacity.txt"])
	subprocess.call(["cp", "./input/planck/opacity_X0000Y0986Z0014.txt", "./LCFiles/kappa_p.txt"])
elif D == 2:
	utils.remesh_evolv_CSM(options.tinj, r_out, CSM_file, file_cc, Ncell=1000)
	i = options.stellar_model.find('Msun')
	Msun = options.stellar_model[i-2:i]
	opac_name = "./input/rosseland/opacity_Icn"+Msun+"Msun.txt"
	subprocess.call(["cp", opac_name, "./LCFiles/opacity.txt"])
	opac_name = "./input/planck/opacity_Icn"+Msun+"Msun.txt"
	subprocess.call(["cp", opac_name, "./LCFiles/kappa_p.txt"])
	

# extract the ejecta parameters
Mej, n, delta, CSM_mass = utils.calculate_ejecta(file_cc, profile_at_cc, CSM_file, D)

# store parameters
with open('params/params.dat', mode='w') as f:
	s = '#The latest parameters used in the calculation are listed.\n'
	f.write(s)
	s = 'Mej = {:.2f} Msun\nn = {:.2f}\nfinj = {:.2f}\ntinj = {:.2f} yr\nmesa model = '.format(Mej, n, options.finj, options.tinj)+options.stellar_model+'\n'+'CSM mass = {:.2f} Msun'.format(CSM_mass)
	f.write(s)

# if multi-band is called, generate frequency-dependent opacity table as well
if options.calc_multiband and D == 0:
	op_freq_dir = 'LCFiles/opacity_frq'
	subprocess.call(["mkdir", "-p", op_freq_dir])
	gen_op_frq.gen_op_frq(Y_He, op_freq_dir)

# calculate light curve
for Eej in options.Eej:
	# luminosity at shock
	dir_Lnu = "LCFiles/SpecFiles_"+str(Eej)
	shock_file = 'LCFiles/shock_output_'+'Mni{:.3f}_'.format(options.Mni)+'tinj{:.2f}_'.format(options.tinj)+str(Eej)+'erg.txt'
#	lightcurve.shock(Eej, Mej*MSUN, options.Mni*MSUN, n, delta, CSM_file, shock_file, D)

	# radiation transfer
	# bolometric light curve
	IIn_lc_file = 'LCFiles/IIn_lightcurve_'+'Mni{:.3f}_'.format(options.Mni)+'tinj{:.2f}_'.format(options.tinj)+str(Eej)+'erg.txt'
	# multi-band light curve if requested
	if options.calc_multiband:
		IIn_lc_band_file = 'LCFiles/IIn_lightcurve_'+'Mni{:.3f}_'.format(options.Mni)+'tinj{:.2f}_'.format(options.tinj)+str(Eej)+'erg_mag.txt'
		subprocess.call(["rm", "-r", dir_Lnu])
		subprocess.call(["mkdir", dir_Lnu])
	else:
		IIn_lc_band_file = ''
	lightcurve.transfer(Eej, Mej*MSUN, options.Mni*MSUN, n, delta, r_out, CSM_file, shock_file, IIn_lc_file, IIn_lc_band_file, dir_Lnu, D)

	# obtain peak luminosity and rise/decay time in days
	# the rise (decay) time is defined by between peak time and the time when the luminosity first rises(decays) to 1% of the peak.
	try:
		utils.extract_peak_and_rise_time(IIn_lc_file, frac=0.01)
	except ValueError:
		pass
