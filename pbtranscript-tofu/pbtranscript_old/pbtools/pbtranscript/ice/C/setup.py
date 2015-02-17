from distutils.core import setup
from distutils.extension import Extension
#from Cython.Distutils import build_ext

ext_modules = [Extension("findECE", ["findECE.c"])]
#ext_modules = [Extension("findECE", ["findECE.pyx"])]

setup(
		name = 'ICE',
        version = '1.1.6',
        description='Isoform-level Clustering for PacBio Transcriptome Data',
        author_email='etseng@pacificbiosciences.com',
        author='Elizabeth Tseng',
        url='http://github.com/PacificBiosciences/cDNA_primer',

        packages = ['Liztools'],
        package_dir = {'Liztools': 'Liztools_pbdagcon/src/Liztools'},
        #cmdclass = {'build_ext': build_ext},
        ext_modules = ext_modules, 
        py_modules = ['basQV','gcon_Liztools',\
            'get_noFL', 'iCEC', 'init_ICE', 'post_ICE', \
            'make_input_fasta_fofn', 'miscBio', \
            'pClique', \
            'findECE', \
            'hmmer_wrapper', \
            'post_quiver_pick_goodones', \
            'quiver_consensus_jchin_way', \
            'include_nonFL_reads', \
            'run_ICE', \
            'run_partial_uc', \
            'run_Quiver_postICE', \
            'utils3'],
        scripts=['hmmer_wrapper.py', 'PBMATRIX.txt', 'include_nonFL_reads.py', 'run_partial_uc.py', 'run_Quiver_postICE.py', 'post_quiver_pick_goodones.py', 'make_input_fasta_fofn.py', 'run_ICE.py', 'split_noFL.sh', 'get_FL_membership.py'],
        data_files=[('example', ['example/c1/g_consensus.fa', 'example/c1/in.fa'])]
)

