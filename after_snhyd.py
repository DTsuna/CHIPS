from __future__ import print_function
from optparse import OptionParser
import sys

# our modules
import utils
import lightcurve
from TOPS import gen_op_tbl 


def parse_command_line():
	parser = OptionParser(
		description = '''Execution script. e.g.,\n
		python run.py --zams-m 15 --zams-z 1 
		'''
	)
	parser.add_option("--zams-m", metavar = "float", type = "float", help = "Initial mass in units of solar mass.")
	parser.add_option("--zams-z", metavar = "float", type = "float", help = "Initial metallicity in units of solar metallicity (Z=0.014).")

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

# extract the ejecta parameters
Mej, n, delta = utils.calculate_ej_from_mesa(file_cc)
Eexp = 1e51

# outer extent of the CSN to feed into the LC calculation
r_out = 9.9e15
CSM_file = 'inp-data/CSM.txt'
Y_He = utils.remesh_CSM(r_out, 'snhydOutput/atCCSN.txt', CSM_file, file_me)

# obtain opacity 
opacity_file = 'inp-data/opacity.txt'
gen_op_tbl.gen_op_tbl(Y_He, opacity_file)

# luminosity at shock
shock_file = 'inp-data/shock_output.txt'
lightcurve.shock(Eexp, Mej*1.99e33, n, delta, CSM_file, shock_file)

# radiation transfer
lc_file = 'outp-data/lightcurve.txt'
lightcurve.transfer(r_out, CSM_file, shock_file, lc_file)

# obtain peak luminosity and rise/decay time in days
# the rise (decay) time is defined by between peak time and the time when the luminosity first rises(decays) to 1% of the peak.
utils.extract_peak_and_rise_time(LC_output, frac=0.01)
