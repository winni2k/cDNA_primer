from distutils.core import setup
from distutils.extension import Extension
#from Cython.Distutils import build_ext

ext_modules = [Extension("BioReaders", ["BioReaders.c"])]

setup(
		name = 'BioReaders',
		ext_modules = ext_modules
)

