from __future__ import print_function
import glob
from optparse import OptionParser
import sys

# our modules
import utils
import lightcurve
from TOPS import gen_op_tbl 


def parse_command_line():
	parser = OptionParser(
		description = '''Execution script. e.g.,\n
		python run.py --zams-m 15 --zams-z 1 --profile-at-cc snhydOutput/intermediate??yr.txt
		'''
	)
	parser.add_option("--zams-m", metavar = "float", type = "float", help = "Initial mass in units of solar mass.")
	parser.add_option("--zams-z", metavar = "float", type = "float", help = "Initial metallicity in units of solar metallicity (Z=0.014).")
	parser.add_option("--profile-at-cc", metavar = "filename", type = "string", help = "The file with the profile at core collapse.")
	parser.add_option("--analytical-CSM", action = "store_true", default=False, help = "Calibrate CSM by analytical profile given in Tsuna et al (2021). The adiabatic CSM profile is extrapolated to the inner region, correcting the profile obtained from adiabatic calculation that includes artificial shock-compression.")

	options, filenames = parser.parse_args()
	available_masses = [13.,14.,15.,16.,17.,18.,19.,20.,22.,24.,26.,28.,30.]
	available_mesa_models = [(mass, 1.) for mass in available_masses]
	if (options.zams_m,options.zams_z) not in available_mesa_models:
		raise ValueError('stellar model not available')

	return options, filenames

# get command line arguments
options, filenames = parse_command_line()


file_cc = 'mesa_models/'+str(int(options.zams_m))+'Msun_Z'+str(0.014*options.zams_z)+'_preccsn.data'
file_me = file_cc 

#################################################################
#								#
#		IIn light curve model of TS20			#
#								#
#################################################################

# outer extent of the CSN to feed into the LC calculation
r_out = 3e16
CSM_file = 'inp-data/CSM.txt'
Y_He, r_edge = utils.remesh_CSM(r_out, options.profile_at_cc, CSM_file, file_me, options.analytical_CSM)

# extract the ejecta parameters
Mej, n, delta = utils.calculate_ejecta(file_cc, options.profile_at_cc, r_edge)
Eexps = [1e51, 3e51, 1e52]

# obtain opacity 
opacity_file = 'inp-data/opacity.txt'
gen_op_tbl.gen_op_tbl(Y_He, opacity_file)

for Eexp in Eexps:
	# luminosity at shock
	shock_file = 'inp-data/shock_output_'+str(Eexp)+'erg.txt'
	lightcurve.shock(Eexp, Mej*1.99e33, n, delta, CSM_file, shock_file)

	# radiation transfer
	IIn_lc_file = 'outp-data/IIn_lightcurve_'+str(Eexp)+'erg.txt'
	lightcurve.transfer(r_out, CSM_file, shock_file, IIn_lc_file)

	# obtain peak luminosity and rise/decay time in days
	# the rise (decay) time is defined by between peak time and the time when the luminosity first rises(decays) to 1% of the peak.
	utils.extract_peak_and_rise_time(IIn_lc_file, frac=0.01)

# obtain light curve at mass eruption
mass_eruption_lc_file = 'outp-data/mass_eruption_lightcurve.txt'
utils.get_mass_eruption_lightcurve(mass_eruption_lc_file)
