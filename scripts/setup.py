from distutils.core import setup
from distutils.extension import Extension
#from Cython.Distutils import build_ext

ext_modules = [Extension("BioReaders", ["BioReaders.c"]), \
        Extension("c_branch", ["c_branch.c"]), \
        Extension("intersection_unique", ["intersection_unique.c"])]

setup(
		name = 'c_cDNA_primer',
		ext_modules = ext_modules,
        author_email='etseng@pacificbiosciences.com',
        author='Elizabeth Tseng'
)

