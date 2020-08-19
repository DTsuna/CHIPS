from __future__ import print_function
import fileinput
from optparse import OptionParser
import subprocess

def parse_command_line():
	parser = OptionParser(
		description = 'execution script'
	)
	parser.add_option("--zams-m", type = "float", help = "Initial mass.")
	parser.add_option("--zams-z", type = "float", help = "Initial metallicity.")
	parser.add_option("--inlist-file", metavar = "filename", help = "Inlist file with the ZAMS mass and metallicity information.")
	options, filenames = parser.parse_args()

	return options, filenames

# get command line arguments
options, filenames = parse_command_line()

# manipulate the file with the given input zams mass and metallicity
for line in fileinput.input(options.inlist_file, inplace=1):
	if 'initial_mass' in line:
		print("      initial_mass = %f" % float(options.zams_m))
	elif 'initial_z' in line:
		print("      initial_z = %f" % float(options.zams_z))
	elif 'Zbase' in line:
		print("      Zbase = %f" % float(options.zams_z))
	else:
		print(line.rstrip())

# compile mesa script
#subprocess.call("./mk")
# run mesa script
#subprocess.call("./rn")

