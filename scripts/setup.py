from distutils.core import setup
from distutils.extension import Extension
#from Cython.Distutils import build_ext
import numpy as np

ext_modules = [Extension("BioReaders", ["BioReaders.c"]), \
        Extension("c_branch", ["c_branch.c"]), \
        Extension("modified_bx_intervals.intersection_unique", ["modified_bx_intervals/intersection_unique.c"])]

setup(
		name = 'c_cDNA_primer',
		ext_modules = ext_modules,
        include_dirs = [np.get_include()],
        author_email='etseng@pacificbiosciences.com',
        author='Elizabeth Tseng'
)

