# Complete History of Interaction-Powered Supernovae (CHIPS)

This is a tool for obtaining the whole history of the progenitors of
interaction-powered transients. Coupled with the MESA stellar evolution
code and several codes implemented by the authors, the user can obtain the
circumstellar matter profile and light curves of the interaction-powered
supernovae, for a selected mass and metallicity of the progenitor star.

## Steps for running:
1. First install MESA
2. git clone this repository in the mesa-rXXXXX
3. Install the light curve script. (The mass eruption script will be compiled at runtime.)
	`python lcsetup.py install`
4. Execute the script run.py with initial mass and metallicity as arguments. For example, to calculate the
evolution of a 15Msun star with metallicity 0.02 (~ solar), given an inlist file containing these parameters run

	`python run.py --zams-m 15 --zams-z 0.02 --inlist-file /path/to/inlist_file`

   There is an example inlist file in the "example_make_pre_ccsn" repository. If using this one, the command is

	`python run.py --zams-m 15 --zams-z 0.02 --inlist-file example_make_pre_ccsn/inlist_common`


References:
1. Kuriyama, Shigeyama (2020), A&A, 635, 127
2. Takei, Shigeyama (2020), PASJ 72, 67
