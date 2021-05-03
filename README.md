# Complete History of Interaction-Powered Supernovae (CHIPS)

CHIPS is a tool for obtaining the whole history of the progenitors of
interaction-powered transients. Coupling the MESA stellar evolution
code and several codes implemented by the authors, the user can obtain the
circumstellar matter profile and light curves of the interaction-powered
supernovae, for a selected mass and metallicity of the progenitor star.

## What can CHIPS do?

CHIPS can generate a realistic CSM from a model-agnostic mass eruption calculation (reference 1), which can serve as a reference for observers to compare with various observations of the CSM. It can also generate bolometric light curves from CSM interaction (reference 2), which can also be compared with observed light curves. The calculation of mass eruption and light curve takes from several hours to half a day on modern CPUs.

## Pre-reqs

The requirement for running CHIPS is quite minimal. One needs only the standard gcc and gfortran, and python3 with numpy installed.

## Steps for running CHIPS:
### With stellar evolution by MESA:
1. First download and install MESA. Please see <http://mesa.sourceforge.net/prereqs.html> for details on how to install.
2. Also install mesa_reader, a module that lets us easily extract data from the mesa output. This is used in some of our Python scripts. Please see <http://mesa.sourceforge.net/output.html> for details.
3. git clone this repository in the mesa-rXXXXX/star repository
4. Copy the test suite from the MESA that you installed. This is to ensure mesa version compatibilities

	`cp -r $MESA_DIR/star/test_suite/example_make_pre_ccsn/src $MESA_DIR/star/CHIPS/example_make_pre_ccsn`

5. Compile the scripts for the light curve.

	`python lcsetup.py install`

6. Execute the script run.py with initial mass and metallicity as arguments. For example, to calculate the history of a star with ZAMS mass 15 Msun and metallicity 1 Zsun (here assumed to be 0.014), with mass eruption occurring 5 years before core-collapse, run

	`python run.py --zams-m 15 --zams-z 1 --delta-t 5 --inlist-file /path/to/inlist_file`

### Skipping the MESA calculation
Alternatively, one can use the MESA pre-SN models generated by us to completely skip the MESA calculation. Our models cover stars of solar metallicity with ZAMS mass range 13-30 Msun, with 1 Msun interval up to 20 Msun and 2 Msun interval from 20 to 30 Msun.

1. The pre-SN models are in a zip file in the directory input/mesa_models/. Once you un-zip this file, you will find MESA data files with the naming showing the mass and metallicity at ZAMS.
2. Clone this repository at any place, and compile the scripts for the light curve.

	`python lcsetup.py install`

3. We also need to install mesa_reader, as it is used in some of our Python scripts. Please see <http://mesa.sourceforge.net/output.html> for details on how to install.
4. Execute run.py with the argument skip-mesa to use our stellar models, e.g.

	`python run.py --zams-m 15 --zams-z 1 --delta-t 5 --skip-mesa`


### Using an already calculated mass CSM calculation
While the mass eruption calculation is running, files with names intermediate??yr.txt, which record the envelope profile ?? years after energy injection, are created. If these files are available, for calculations of the light curve one can skip the mass-eruption calculation and obtain the light curve using the code after_snhyd.py.

	`python3 after_snhyd.py --zams-m 15 --zams-z 1 --profile-at-cc EruptionFiles/intermediate10yr.txt`

The profile-at-cc argument sets the interval between mass eruption and core-collapse. The files are produced at an interval of 1 year.

### Using an analytical CSM model
We strongly advise to use an analytical CSM model (reference 3) that corrects the artifical shock-compressions that arise from the adiabatic mass eruption code. This can be done with adding an argument --analytical-CSM.


## References:
1. Kuriyama, Shigeyama (2020), A&A, 635, 127 (for mass eruption)
2. Takei, Shigeyama (2020), PASJ 72, 67 (for light curve)
3. Tsuna, Takei, Kuriyama, Shigeyama (2021), arXiv: 2104.03694
