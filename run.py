import fileinput
from optparse import OptionParser
import subprocess
import sys

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

# edit the file with given input zams and metallicity
for line in fileinput.input(options.inlist_file, inplace=1):
	if 'initial_mass' in line:
		print "      initial_mass = %f\n" % options.zams_m,
	elif 'initial_z' in line:
		print "      initial_z = %f\n" % options.zams_z,
	elif 'Zbase' in line:
		print "      Zbase = %f\n" % options.zams_z,
	else:
		print line,

# compile mesa script
#subprocess.call("./mk")
# run mesa script
#subprocess.call("./rn")


