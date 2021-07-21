from __future__ import print_function
import glob
from optparse import OptionParser
import sys

# our modules
from utils import utils
from utils import mag
import lightcurve
from input.TOPS import gen_op_tbl 


def parse_command_line():
	parser = OptionParser(
		description = '''Execution script. e.g.,\n
		python run.py --zams-m 15 --zams-z 1 --profile-at-cc EruptionFiles/intermediate??yr.txt
		'''
	)
	parser.add_option("--zams-m", metavar = "float", type = "float", help = "Initial mass in units of solar mass.")
	parser.add_option("--zams-z", metavar = "float", type = "float", help = "Initial metallicity in units of solar metallicity (Z=0.014).")
	parser.add_option("--profile-at-cc", metavar = "filename", type = "string", help = "The file with the profile at core collapse.")
	parser.add_option("--analytical-CSM", action = "store_true", default=False, help = "Calibrate CSM by analytical profile given in Tsuna et al (2021). The adiabatic CSM profile is extrapolated to the inner region, correcting the profile obtained from adiabatic calculation that includes artificial shock-compression.")
	parser.add_option("--steady-wind", metavar = "string", type = "string", default='RSGwind', help = "Set how the steady wind CSM is attached to the erupted material. Must be 'attach' or 'RSGwind'.")

	options, filenames = parser.parse_args()
	available_masses = [13.,14.,15.,16.,17.,18.,19.,20.,22.,24.,26.,28.,30.]
	available_mesa_models = [(mass, 1.) for mass in available_masses]
	if (options.zams_m,options.zams_z) not in available_mesa_models:
		raise ValueError('stellar model not available')
	if options.profile_at_cc == None:
		raise ValueError('the density profile at core-collapse (EruptionFiles/intermediate??.txt) is a required argument')

	return options, filenames

# get command line arguments
options, filenames = parse_command_line()


file_cc = 'input/mesa_models/'+str(int(options.zams_m))+'Msun_Z'+str(0.014*options.zams_z)+'_preccsn.data'
file_me = file_cc 

#################################################################
#								#
#		IIn light curve model of TS20			#
#								#
#################################################################

# outer extent of the CSN to feed into the LC calculation
r_out = 3e16
CSM_file = 'LCFiles/CSM.txt'
Y_He, r_edge = utils.remesh_CSM(r_out, options.profile_at_cc, CSM_file, file_me, analytical_CSM = options.analytical_CSM, steady_wind = options.steady_wind)

# extract the ejecta parameters
Mej, n, delta = utils.calculate_ejecta(file_cc, options.profile_at_cc, CSM_file, r_edge)
Eexps = [1e51, 3e51, 1e52]

# obtain opacity 
opacity_file = 'LCFiles/opacity.txt'
gen_op_tbl.gen_op_tbl_sct(Y_He, opacity_file)
opacity_file = 'LCFiles/kappa_p.txt'
gen_op_tbl.gen_op_tbl_abs(Y_He, opacity_file)

for Eexp in Eexps:
	# luminosity at shock
	shock_file = 'LCFiles/shock_output_'+str(Eexp)+'erg.txt'
	lightcurve.shock(Eexp, Mej*1.99e33, n, delta, CSM_file, shock_file)

	# radiation transfer
	IIn_lc_file = 'LCFiles/IIn_lightcurve_'+str(Eexp)+'erg.txt'
	lightcurve.transfer(Eexp, Mej*1.99e33, n, delta, r_out, CSM_file, shock_file, IIn_lc_file)

	# obtain peak luminosity and rise/decay time in days
	# the rise (decay) time is defined by between peak time and the time when the luminosity first rises(decays) to 1% of the peak.
	utils.extract_peak_and_rise_time(IIn_lc_file, frac=0.01)
	mag.calc_mag(IIn_lc_file)

# obtain light curve at mass eruption
mass_eruption_lc_file = 'LCFiles/mass_eruption_lightcurve.txt'
utils.get_mass_eruption_lightcurve(mass_eruption_lc_file)
