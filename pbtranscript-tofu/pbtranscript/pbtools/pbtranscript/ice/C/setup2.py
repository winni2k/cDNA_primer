from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("c_basQV", ["c_basQV.pyx"], language="c++"), \
#        Extension("c_Prob", ["ProbHandler.pyx"], language="c++"), \
        Extension("ProbModel", ["ProbModel.pyx"], language="c++"),\
        Extension("findECE", ["findECE.pyx"])]

setup(
		name = 'c_ICE',
		cmdclass = {'build_ext': build_ext},
		ext_modules = ext_modules
)

