from optparse import OptionParser
import subprocess

# our modules
from utils import utils
from input.TOPS import gen_op_tbl 
from input.TOPS_multigroup import gen_op_frq
import lightcurve


def parse_command_line():
	parser = OptionParser(
		description = '''Execution script. e.g.,\n
		python run.py --Eej 1e51 --stellar-model input/mesa_models/15Msun_Z0.014_preccsn.data --profile-at-cc EruptionFiles/intermediate??yr.txt --analytical-CSM
		'''
	)
	parser.add_option("--Eej", metavar = "float", type = "float", action = "append", help = "Explosion energy in erg. This option can be given multiple times (default: 1e51, 3e51, 1e52).")
	parser.add_option("--Mni", metavar = "float", type = "float", default = 0., help = "Radioactive nickel 56 mass, in units of solar mass (default: 0 Msun).")
	parser.add_option("--stellar-model", metavar = "filename", help = "Path to the input stellar model (required). This should be one of the stellar model files created after running MESA (which usually end with '.data'.). If --run-mesa is called, this needs to be the stellar model file that you want to provide as input of the CHIPS code (e.g. the file provided by the input 'filename_for_profile_when_terminate' in one of the inlist files.).")
	parser.add_option("--profile-at-cc", metavar = "filename", type = "string", help = "The file with the profile at core collapse.")
	parser.add_option("--analytical-CSM", action = "store_true", default=False, help = "Calibrate CSM by analytical profile given in Tsuna et al (2021). The adiabatic CSM profile is extrapolated to the inner region, correcting the profile obtained from adiabatic calculation that includes artificial shock-compression.")
	parser.add_option("--steady-wind", metavar = "string", type = "string", default='RSGwind', help = "Specify how the steady wind CSM is attached to the erupted material. Must be 'attach' or 'RSGwind' (default: RSGwind). 'attach' simply connects a wind profile to the outermost cell profile, while 'RSGwind' smoothly connects a red supergiant wind to the erupted material.")
	parser.add_option("--calc-multiband", action = "store_true", default=False, help = "Additionally conduct ray-tracing calculations to obtain multi-band light curves (default: false). This calculation is computationally heavier than obtaining just the bolometric light curve.")

	options, filenames = parser.parse_args()
	# sanity checks
	assert options.stellar_model is not None, "A stellar model file for input to CHIPS has to be provided."
	if options.profile_at_cc is None:
		raise ValueError('the density profile at core-collapse (EruptionFiles/intermediate??.txt) is a required argument')
	# set default value if explosion energy is empty
	if not options.Eej:
		options.Eej = [1e51, 3e51, 1e52]

	return options, filenames

# get command line arguments
options, filenames = parse_command_line()

file_cc = options.stellar_model 

# get tinj from profile_at_cc name
tinj = float(options.profile_at_cc.split('yr.txt')[0][-2])

#################################################################
#								#
#		IIn light curve model of TS20			#
#								#
#################################################################

# discriminate progenitor type
# D=0: stars with hydrogen-rich envelope (SNType='IIn')
# D=1: stripped stars with helium-rich outer layer (SNType='Ibn')
# D=2: stripped stars with helium-poor outer layer (SNType='Icn')
D, SNType = utils.discriminantProgModel(file_cc)

# store abundance
utils.genAbundanceTable(file_cc)

# generate tables for mu, T.
lightcurve.opacTable(D)


# outer extent of the CSN to feed into the LC calculation
r_out = 3e16
# remesh CSM in order to correct for shocks in the hydro simulation and extend to r_out.
CSM_file = 'LCFiles/CSM.txt'
Y_He = utils.remesh_CSM(r_out, options.profile_at_cc, CSM_file, file_cc, analytical_CSM = options.analytical_CSM, steady_wind = options.steady_wind)

# extract the ejecta parameters
Mej, n, delta, CSM_mass = utils.calculate_ejecta(file_cc, options.profile_at_cc, CSM_file, D)

with open('params/params_after_eruption.dat', mode='w') as f:
	s = '#The latest parameters used in the calculation are listed.\n'
	f.write(s)
	s = 'Mej = {:.2f} Msun\nn = {:.2f}\nmesa model = '.format(Mej, n)+options.stellar_model+'\n'+'profile at core collapse = '+options.profile_at_cc+'\n'+'CSM mass = {:.2f} Msun'.format(CSM_mass)
	f.write(s)

# only IIn is supported
assert SNType == 'IIn'

# obtain opacity 
opacity_file = 'LCFiles/opacity.txt'
gen_op_tbl.gen_op_tbl_sct(Y_He, opacity_file)
opacity_file = 'LCFiles/kappa_p.txt'
gen_op_tbl.gen_op_tbl_abs(Y_He, opacity_file)
# if multi-band is called, generate frequency-dependent opacity table as well
if options.calc_multiband:
	op_freq_dir = 'LCFiles/opacity_frq'
	subprocess.call(["mkdir", "-p", op_freq_dir])
	gen_op_frq.gen_op_frq(Y_He, op_freq_dir)

# calculate light curve
for Eej in options.Eej:
	# luminosity at shock
	shock_file = 'LCFiles/{}_shock_output_'.format(SNType)+'Mni{:.3f}_'.format(options.Mni)+'tinj{:.2f}_'.format(tinj)+str(Eej)+'erg.txt'
	lightcurve.shock(Eej, Mej*1.99e33, options.Mni, n, delta, CSM_file, shock_file, D)

	# radiation transfer
	# bolometric light curve
	IIn_lc_file = 'LCFiles/{}_lightcurve_'.format(SNType)+'Mni{:.3f}_'.format(options.Mni)+'tinj{:.2f}_'.format(tinj)+str(Eej)+'erg.txt'
	dir_Lnu = "LCFiles/SpecFiles_"+str(Eej)
	# multi-band light curve if requested
	if options.calc_multiband:
		IIn_lc_band_file = 'LCFiles/{}_lightcurve_'.format(SNType)+'Mni{:.3f}_'.format(options.Mni)+'tinj{:.2f}_'.format(tinj)+str(Eej)+'erg_mag.txt'
		subprocess.call(["rm", "-r", dir_Lnu])
		subprocess.call(["mkdir", dir_Lnu])
	else:
		IIn_lc_band_file = ''
	lightcurve.transfer(Eej, Mej*1.99e33, options.Mni, n, delta, r_out, CSM_file, shock_file, IIn_lc_file, IIn_lc_band_file, dir_Lnu, D)

	# obtain peak luminosity and rise/decay time in days
	# the rise (decay) time is defined by between peak time and the time when the luminosity first rises(decays) to 1% of the peak.
	try:
		utils.extract_peak_and_rise_time(IIn_lc_file, frac=0.01)
	except ValueError:
		pass
