from distutils.core import setup, Extension

module1 = Extension('lightcurve',
                    sources = ['lcwrapper.c', 'src/shock.c', 'transfer/src/transfer.c'],
                    )

setup(name = 'lightcurve', version = '1.0.0', ext_modules = [module1])
