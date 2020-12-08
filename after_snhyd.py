from __future__ import print_function
import glob
import sys

# our modules
import utils
import lightcurve
from TOPS import gen_op_tbl 

lctest_dir = 'inp-data/' 

# find data file at mass eruption and core collapse. 
# FIXME we set the mass eruption to 5 years before collapse
file_cc = lctest_dir+'pre_ccsn.data'
file_me = lctest_dir+'profile14.data'

#################################################################
#								#
#		IIn light curve model of TS20			#
#								#
#################################################################

# extract the ejecta parameters
print("extracting ejecta parameters...", file=sys.stderr)
Mej, n, delta = utils.calculate_ej_from_mesa(file_cc)
Eexp = 1e51

# outer extent of the CSN to feed into the LC calculation
print("remeshing CSM...", file=sys.stderr)
r_out = 9.9e15
CSM_file = 'inp-data/CSM.txt'
Y_He = utils.remesh_CSM(r_out, 'snhydOutput/atCCSN.txt', CSM_file, file_me)

# obtain opacity 
print("obtaining opacity..,", file=sys.stderr)
opacity_file = 'inp-data/opacity.txt'
gen_op_tbl.gen_op_tbl(Y_He, opacity_file)

# luminosity at shock
print("calculating light curve...", file=sys.stderr)
shock_file = 'inp-data/shock_output.txt'
lightcurve.shock(Eexp, Mej*1.99e33, n, delta, CSM_file, shock_file)

# radiation transfer
lc_file = 'outp-data/lightcurve.txt'
lightcurve.transfer(r_out, CSM_file, shock_file, lc_file)

# obtain peak luminosity and rise/decay time in days
# the rise (decay) time is defined by between peak time and the time when the luminosity first rises(decays) to 1% of the peak.
peakL, riset, decayt = utils.extract_peak_and_rise_time(LC_output, frac=0.01)
