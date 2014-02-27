import distutils
from distutils.core import setup, Extension
import distutils.ccompiler
import glob
import os
import os.path
import sys
import numpy as np

if sys.argv[1] == "clean":
    print "Cleaning files..."
    os.system("rm -rf ./build/")
    os.system("rm -rf misopy/pyx/*.c")
    os.system("rm -rf misopy/pyx/*.so")

# Extract long description of MISO from README
long_description = open('README').read()

if sys.version_info > (3, 0):
    options["use_2to3"] = True

# This forces distutils to place the data files
# in the directory where the Py packages are installed
# (usually 'site-packages'). This is unfortunately
# required since there's no good way to retrieve
# data_files= from setup() in a platform independent
# way.
from distutils.command.install import INSTALL_SCHEMES
for scheme in INSTALL_SCHEMES.values():
        scheme['data'] = scheme['purelib']

##
## Definition of the current version
##
MISO_VERSION = "0.6.0"

##
## Generate a __version__.py attribute
## for the module
##
with open("./misopy/__init__.py", "w") as version_out:
      version_out.write("__version__ = \"%s\"\n" %(MISO_VERSION))

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
      
# Include our headers and numpy's headers
include_dirs = [os.path.join(CURRENT_DIR, "include")]

##
## Extension modules
##
# Single-end scoring functions
single_end_ext = Extension("misopy.pyx.miso_scores_single",
                           ["misopy/pyx/miso_scores_single.pyx"],
                           libraries=["m"])
#                           libraries=["m"],
#                           include_dirs=include_dirs)

# Paired-end scoring functions
paired_end_ext = Extension("misopy.pyx.miso_scores_paired",
                           ["misopy/pyx/miso_scores_paired.pyx"],
                           libraries=["m"])
#                           libraries=["m"],                           
#                           include_dirs=include_dirs)

# Add sampler routine here...
# ....


# Statistics functions
stat_helpers_ext = Extension("misopy.pyx.stat_helpers",
                             ["misopy/pyx/stat_helpers.pyx"],
                             libraries=["m"])

# Matrix functions
matrix_utils_ext = Extension("misopy.pyx.matrix_utils",
                             ["misopy/pyx/matrix_utils.pyx"],
                             libraries=["m"])


#                             libraries=["m"],
#                             include_dirs=include_dirs)
# Lapack functions extension
cc = distutils.ccompiler.new_compiler()
defines = []
if cc.has_function('rintf(1.0);rand', includes=['math.h', 'stdlib.h'],
                   libraries=['m']):
    defines.append(('HAVE_RINTF', '1'))
if cc.has_function('finite(1.0);rand', includes=['math.h', 'stdlib.h']):
    defines.append(('HAVE_FINITE', '1'))
if cc.has_function('expm1(1.0);rand', includes=['math.h', 'stdlib.h'],
                   libraries=['m']):
    defines.append(('HAVE_EXPM1', '1'))
if cc.has_function('rint(1.0);rand', includes=['math.h', 'stdlib.h'], 
                   libraries=['m']):
    defines.append(('HAVE_RINT', '1'))
if cc.has_function('double log2(double) ; log2(1.0);rand', 
                   includes=['math.h', 'stdlib.h'], libraries=['m']):
    defines.append(('HAVE_LOG2', '1'))
if cc.has_function('logbl(1.0);rand', includes=['math.h', 'stdlib.h'],
                   libraries=['m']):
    defines.append(('HAVE_LOGBL', '1'))
if cc.has_function('snprintf(0, 0, "");rand', 
                   includes=['stdio.h', 'stdlib.h']):
    defines.append(('HAVE_SNPRINTF', '1'))
if cc.has_function('log1p(1.0);rand', includes=['math.h', 'stdlib.h'],
                   libraries=['m']):
    defines.append(('HAVE_LOG1P', '1'))
if cc.has_function('double round(double) ; round(1.0);rand', 
                   includes=['math.h', 'stdlib.h'], libraries=['m']):
    defines.append(('HAVE_ROUND', '1'))
if cc.has_function('double fmin(double, double); fmin(1.0,0.0);rand', 
                   includes=['math.h', 'stdlib.h'], libraries=['m']):
    defines.append(('HAVE_FMIN', '1'))
# Source files
c_source_dir = os.path.join(CURRENT_DIR, "src")
lapack_sources = \
  glob.glob(os.path.join(c_source_dir, "lapack", "*.c"))
f2c_sources = \
  glob.glob(os.path.join(c_source_dir, "f2c", "*.c"))
blas_sources = \
  glob.glob(os.path.join(c_source_dir, "blas", "*.c"))
# Include numpy headers
#all_c_sources = \
#  lapack_sources + blas_sources + f2c_sources
#all_c_sources 

#lapack_ext = Extension("misopy.pyx.lapack",
#                       all_c_sources + ["misopy/pyx/lapack.pyx"],
#                       libraries=["m"],                       
#                       include_dirs=include_dirs)
#                       define_macros=defines)

# pyx/c extensions to MISO
miso_extensions = [single_end_ext,
                   paired_end_ext,
                   stat_helpers_ext,
                   matrix_utils_ext]

#                   lapack_ext]

##
## Handle creation of source distribution. Here we definitely
## need to use Cython. This creates the *.c files from *.pyx
## files.
##
cmdclass = {}
from distutils.command.sdist import sdist as _sdist
class sdist(_sdist):
    """
    Override sdist command to use cython
    """
    def run(self):
        try:
            from Cython.Build import cythonize
        except ImportError:
            raise Exception, "Cannot create source distribution without Cython."
        print "Cythonizing"
        extensions = cythonize(miso_extensions)
        _sdist.run(self)
cmdclass['sdist'] = sdist


##
## Handle Cython sources. Determine whether or not to use Cython
##
extensions = []
def no_cythonize(extensions, **_ignore):
    new_extensions = []
    for extension in extensions:
        # Make copy of extension so we do't modify
        # it destructively
        ext_copy = \
          Extension(extension.name,
                    extension.sources,
                    include_dirs=extension.include_dirs,
                    library_dirs=extension.library_dirs)
        sources = []
        for sfile in ext_copy.sources:
            path, ext = os.path.splitext(sfile)
            new_sfile = sfile
            if ext in ('.pyx', '.py'):
                if extension.language == 'c++':
                    ext = '.cpp'
                else:
                    ext = '.c'
                new_sfile = path + ext
            sources.append(new_sfile)
        ext_copy.sources[:] = sources
        new_extensions.append(ext_copy)
    return new_extensions


# Whether or not to use Cython
USE_CYTHON = False

try:
    from Cython.Distutils import build_ext
    from Cython.Build import cythonize
    USE_CYTHON = True
except ImportError:
    USE_CYTHON = False


## Force to not use Cython, unless we're making
## a source distribution with 'sdist'
if sys.argv[1] == "sdist":
    USE_CYTHON = True
else:
    USE_CYTHON = False

if USE_CYTHON:
    print "Using Cython."
    extensions = cythonize(miso_extensions)
    cmdclass.update({'build_ext': build_ext})
else:
    extensions = no_cythonize(miso_extensions)
    from distutils.command import build_ext
    print "Not using Cython."


setup(name = 'misopy',
      ##
      ## CRITICAL: When changing version, remember
      ## to change version in __version__.py
      ##
      version = MISO_VERSION,
      description = "Mixture of Isoforms model (MISO) " \
                    "for isoform quantitation using RNA-Seq",
      long_description = long_description,
      author = 'Yarden Katz,Gabor Csardi',
      author_email = 'yarden@mit.edu,gcsardi@stat.harvard.edu',
      maintainer = 'Yarden Katz',
      maintainer_email = 'yarden@mit.edu',
      url = 'http://genes.mit.edu/burgelab/miso/',
      cmdclass = cmdclass,
      ext_modules = extensions,
      include_dirs = [np.get_include()],
      # Tell distutils to look for pysplicing in the right directory
      package_dir = {'pysplicing': 'pysplicing/pysplicing'},
      packages = ['misopy',
                  'misopy.pyx',
                  'misopy.sashimi_plot',
                  'misopy.sashimi_plot.plot_utils',
                  'misopy.internal_tests',
                  'pysplicing'],
      entry_points={
          'console_scripts': [
              'module_availability = misopy.module_availability:main',
              'sam_to_bam = misopy.sam_to_bam:main',
              'run_events_analysis.py = misopy.run_events_analysis:main',
              'run_miso.py = misopy.run_miso:main',
              'exon_utils = misopy.exon_utils:main',
              'pe_utils = misopy.pe_utils:main',
              'filter_events = misopy.filter_events:main',
              'test_miso = misopy.test_miso:main',
              'miso_zip = misopy.miso_zip:main',
              'miso = misopy.miso:main',
              'compare_miso = misopy.compare_miso:main',
              'summarize_miso = misopy.summarize_miso:main',
              'index_gff = misopy.index_gff:main',
              'miso_pack = misopy.miso_pack:main',
              # sashimi_plot scripts
              'plot.py = misopy.sashimi_plot.plot:main',
              'sashimi_plot = misopy.sashimi_plot.sashimi_plot:main'
              ],
          },                 
      data_files = [('misopy/settings',
                     ['misopy/settings/miso_settings.txt',
                      'misopy/sashimi_plot/settings/sashimi_plot_settings.txt']),
                    ('misopy/test-data/sam-data',
                     ['misopy/test-data/sam-data/c2c12.Atp2b1.sam']),
                    ('misopy/sashimi_plot/test-data', 
                      ['misopy/sashimi_plot/test-data/events.gff']),
                    ('misopy/sashimi_plot/test-data/miso-data',
                     ['misopy/sashimi_plot/test-data/miso-data/heartKOa/chr17/chr17:45816186:45816265:-@chr17:45815912:45815950:-@chr17:45814875:45814965:-.miso',
                      'misopy/sashimi_plot/test-data/miso-data/heartKOb/chr17/chr17:45816186:45816265:-@chr17:45815912:45815950:-@chr17:45814875:45814965:-.miso',
                      'misopy/sashimi_plot/test-data/miso-data/heartWT1/chr17/chr17:45816186:45816265:-@chr17:45815912:45815950:-@chr17:45814875:45814965:-.miso',
                      'misopy/sashimi_plot/test-data/miso-data/heartWT2/chr17/chr17:45816186:45816265:-@chr17:45815912:45815950:-@chr17:45814875:45814965:-.miso']),
                    ('misopy/sashimi_plot/test-data/bam-data',
                     ['misopy/sashimi_plot/test-data/bam-data/heartKOa.bam.bai',
                      'misopy/sashimi_plot/test-data/bam-data/heartKOa.sorted.bam',
                      'misopy/sashimi_plot/test-data/bam-data/heartKOa.sorted.bam.bai',
                      'misopy/sashimi_plot/test-data/bam-data/heartKOb.sorted.bam',
                      'misopy/sashimi_plot/test-data/bam-data/heartKOb.sorted.bam.bai',
                      'misopy/sashimi_plot/test-data/bam-data/heartWT1.sorted.bam',
                      'misopy/sashimi_plot/test-data/bam-data/heartWT1.sorted.bam.bai',
                      'misopy/sashimi_plot/test-data/bam-data/heartWT2.sorted.bam',
                      'misopy/sashimi_plot/test-data/bam-data/heartWT2.sorted.bam.bai'])],
      # Required modules
      install_requires = [
          "matplotlib",
          "numpy >= 1.5.0",
          "scipy >= 0.9.0",
          "pysam >= 0.6.0"
          ],
      platforms = 'ALL',
      keywords = ['bioinformatics', 'sequence analysis',
                  'alternative splicing', 'RNA-Seq',
                  'probabilistic models', 'bayesian'],
      classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Programming Language :: C',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ]
      )

