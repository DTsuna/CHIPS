from distutils.core import setup, Extension
import glob

shock_src = glob.glob('./shock/src/*')
transfer_src = glob.glob('./transfer/src/*')

module1 = Extension('lightcurve',
                    define_macros = [('NSIZE', '1000'), ('EADD', None)],
                    sources = ['lcwrapper.c']+shock_src+transfer_src,
		    include_dirs = ['./shock/include', './transfer/include'],
                    )

setup(name = 'lightcurve', version = '1.0.0', ext_modules = [module1])
