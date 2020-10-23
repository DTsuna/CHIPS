from distutils.core import setup, Extension
import glob

string = glob.glob('./src/*')

module1 = Extension('lightcurve',
                    sources = ['lcwrapper.c']+string,
		    include_dirs = ['./include'],
                    )

setup(name = 'lightcurve', version = '1.0.0', ext_modules = [module1])
