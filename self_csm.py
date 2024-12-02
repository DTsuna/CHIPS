import glob
from optparse import OptionParser
import os
import sys
import subprocess
import datetime

# our modules
from utils import convert
from utils import utils
from input.TOPS import gen_op_tbl 
from input.TOPS_multigroup import gen_op_frq
import lightcurve


def parse_command_line():
	parser = OptionParser(
		description = '''Execution script. The MESA calculation will be conducted by running the code with ZAMS mass/Z and inlist file as arguments, For example, to calculate the interaction-powered supernova of a ZAMS 15Msun, solar metallicity star with tinj=10yrs, finj=0.5, and explosion energy 1e51 ergs, run the following:\n
		python run.py --tinj 10 --finj 0.5 --Eej 1e51 --stellar-model input/mesa_models/15Msun_Z0.014_preccsn.data --analytical-CSM
		'''
	)
	parser.add_option("--Eej", metavar = "float", type = "float", action = "append", help = "Explosion energy in erg. This option can be given multiple times (default: 1e51, 3e51, 1e52).")
	parser.add_option("--Mni", metavar = "float", type = "float", default = 0., help = "Radioactive nickel 56 mass, in units of solar mass (default: 0 Msun).")
	parser.add_option("--exponent", metavar = "float", type = "float", default = -2., help = "power-law index of CSM (default: stellar wind: -2).")
	parser.add_option("--CSM-mass", metavar = "float", type = "float", help = "CSM mass, in units of solar mass.")
	parser.add_option("--break-point", metavar = "float", type = "float", help = "CSM mass, in units of solar mass.")
	parser.add_option("--stellar-model", metavar = "filename", help = "Path to the input stellar model (required). This should be one of the stellar model files created after running MESA (which usually end with '.data'.). If --run-mesa is called, this needs to be the stellar model file that you want to provide as input of the CHIPS code (e.g. the file provided by the input 'filename_for_profile_when_terminate' in one of the inlist files.).")
	parser.add_option("--run-mesa", action = "store_true", help = "Call to run MESA in this script and get a new stellar model.")
	parser.add_option("--mesa-path", metavar = "string", type = "string", help = "Path to the execution files of MESA.")
	parser.add_option("--steady-wind", metavar = "string", type = "string", default='RSGwind', help = "Specify how the steady wind CSM is attached to the erupted material. Must be 'attach' or 'RSGwind' (default: RSGwind). 'attach' simply connects a wind profile to the outermost cell profile, while 'RSGwind' smoothly connects a red supergiant wind to the erupted material.")
	parser.add_option("--calc-multiband", action = "store_true", default=False, help = "Additionally conduct ray-tracing calculations to obtain multi-band light curves (default: false). This calculation is computationally heavier than obtaining just the bolometric light curve. For now this feature is only enabled for Type IIn cases.")
	parser.add_option("--opacity-table",  metavar = "filename", help = "A custom opacity table used for the mass eruption calculation. If not called, analytical opacity formula (Kuriyama & Shigeyama 20) is used. The opacity file should have the format like the files in the input/rosseland directory.")

	options, filenames = parser.parse_args()

	# sanity checks
	assert options.stellar_model is not None, "A stellar model file for input to CHIPS has to be provided."
	if options.run_mesa:
		assert options.mesa_path is not None, "A valid existing directory has to be given for --mesa-path if --run-mesa is called."
	# set default value if explosion energy is empty
	if not options.Eej:
		options.Eej = [1e51, 3e51, 1e52]
	if options.CSM_mass is None or options.break_point is None:
		raise ValueError("Parameters CSM-mass and break-point need to be provided.")

	return options, filenames


# get command line arguments
options, filenames = parse_command_line()

# store parameters
dt_now = datetime.datetime.now()
s='{:02d}{:02d}{:02d}_{:02d}{:02d}{:02d}'.format(dt_now.year%100,dt_now.month,dt_now.day,dt_now.hour,dt_now.minute,dt_now.second)
with open('params/params_'+s+'.dat', mode='w') as f:
	s = '#The latest parameters used in the calculation are listed.\n'
	f.write(s)
	s = 'power-law index = {:.2f}\nCSM mass = {:.2f} Msun\nNickel mass = {:.3f} Msun\nmesa model = '.format(options.exponent, options.CSM_mass, options.Mni)+options.stellar_model+'\n'
	Eej_str = ['%g erg' % E for E in options.Eej]
	s = s+'Eej = '+str(Eej_str)+'\n'
	s = s+'multi_band = '+str(options.calc_multiband)
	f.write(s)


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
# D=0: stars with hydrogen-rich envelope (SNType='IIn')
# D=1: stripped stars with helium-rich outer layer (SNType='Ibn')
# D=2: stripped stars with helium-poor outer layer (SNType='Icn')
D, SNType = utils.discriminantProgModel(file_cc)

# store abundance
utils.genAbundanceTable(file_cc)

# generate tables for mu, T.
lightcurve.opacTable(D)

#################################################################
#								#
#		Light curve model of TS20			#
#								#
#################################################################


# outer extent of the CSN to feed into the LC calculation
r_out = 3e16
# remesh CSM in order to correct for shocks in the hydro simulation and extend to r_out.
CSM_file = 'LCFiles/CSM.txt'
profile_at_cc = 'EruptionFiles/atCCSN.txt'
profile_at_cc = [options.exponent, options.CSM_mass, options.break_point]
if SNType == 'IIn':
# obtain opacity
	Y_He = utils.remesh_CSM(r_out, profile_at_cc, CSM_file, file_cc, steady_wind=options.steady_wind)
	opacity_file = 'LCFiles/opacity.txt'
	gen_op_tbl.gen_op_tbl_sct(Y_He, opacity_file)
	opacity_file = 'LCFiles/kappa_p.txt'
	gen_op_tbl.gen_op_tbl_abs(Y_He, opacity_file)
elif SNType == 'Ibn':
	utils.remesh_evolv_CSM(options.tinj, r_out, CSM_file, file_cc, Ncell=1000)
	subprocess.call(["cp", "./input/rosseland/opacity_X0000Y0986Z0014.txt", "./LCFiles/opacity.txt"])
	subprocess.call(["cp", "./input/planck/opacity_X0000Y0986Z0014.txt", "./LCFiles/kappa_p.txt"])
elif SNType == 'Icn':
	utils.remesh_evolv_CSM(options.tinj, r_out, CSM_file, file_cc, Ncell=1000)
	i = options.stellar_model.find('Msun')
	Msun = options.stellar_model[i-2:i]
	opac_name = "./input/rosseland/opacity_Icn"+Msun+"Msun.txt"
	subprocess.call(["cp", opac_name, "./LCFiles/opacity.txt"])
	opac_name = "./input/planck/opacity_Icn"+Msun+"Msun.txt"
	subprocess.call(["cp", opac_name, "./LCFiles/kappa_p.txt"])
else:
	raise ValueError('invalid SNType: %s' % SNType)


# extract the ejecta parameters
Mej, n, delta, CSM_mass = utils.calculate_ejecta(file_cc, profile_at_cc, CSM_file, D)
utils.interpolate_self_similar_solution(n, options.exponent)

# if multi-band is called, generate frequency-dependent opacity table as well
if options.calc_multiband and SNType=='IIn':
	op_freq_dir = 'LCFiles/opacity_frq'
	subprocess.call(["mkdir", "-p", op_freq_dir])
	gen_op_frq.gen_op_frq(Y_He, op_freq_dir)

# calculate light curve
for Eej in options.Eej:
	print('Currently working on Eej=%g...' % Eej)

	# luminosity at shock
	dir_Lnu = "LCFiles/SpecFiles_"+str(Eej)
	shock_file = 'LCFiles/{}_shock_output_'.format(SNType)+'Mni{:.3f}_'.format(options.Mni)+'CSM_mass{:.2f}Msun_exponent{:.2f}_breakpoint{:.1e}cm_'.format(options.CSM_mass, options.exponent, options.break_point)+str(Eej)+'erg.txt'
	lightcurve.shock(Eej, Mej, options.Mni, n, delta, -options.exponent, CSM_file, shock_file, D)

	# radiation transfer
	# bolometric light curve
	lc_file = 'LCFiles/{}_lightcurve_'.format(SNType)+'Mni{:.3f}_'.format(options.Mni)+'CSM_mass{:.2f}Msun_exponent{:.2f}_breakpoint{:.1e}cm_'.format(options.CSM_mass, options.exponent, options.break_point)+str(Eej)+'erg.txt'
	# multi-band light curve if requested
	lc_band_file = ''
	if options.calc_multiband:
		if SNType=='IIn':
			lc_band_file = 'LCFiles/{}_lightcurve_'.format(SNType)+'Mni{:.3f}_'.format(options.Mni)+'CSM_mass{:.2f}Msun_exponent{:.2f}_breakpoint{:.1e}cm_'.format(options.CSM_mass, options.exponent, options.break_point)+str(Eej)+'erg_mag.txt'
			subprocess.call(["rm", "-r", dir_Lnu])
			subprocess.call(["mkdir", dir_Lnu])
		else:
			print('Multi-band light curves are currently enabled only for IIn. Skipping for Ibn/Icn...')

	lightcurve.transfer(Eej, Mej, options.Mni, n, delta, r_out, -options.exponent, CSM_file, shock_file, lc_file, lc_band_file, dir_Lnu, D)

	# obtain peak luminosity and rise/decay time in days
	# the rise (decay) time is defined by between peak time and the time when the luminosity first rises(decays) to 1% of the peak.
	try:
		utils.extract_peak_and_rise_time(lc_file, frac=0.01)
	except ValueError:
		pass
