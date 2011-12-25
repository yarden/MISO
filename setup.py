from distutils.core import setup, Extension
import distutils.ccompiler
import glob
import os
import sys

## Test for functions, with a hack to suppress compiler warnings.
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

# new -- prefix for pysplicing
pysplicing_dir = 'pysplicing'

# new
splicingsources = glob.glob(os.path.join(pysplicing_dir, 'src', '*.c'))
lapacksources = glob.glob(os.path.join(pysplicing_dir, 'src', 'lapack', '*.c'))
f2csources = glob.glob(os.path.join(pysplicing_dir, 'src', 'f2c', '*.c'))

#old
#splicingsources = glob.glob(os.path.join('src', '*.c'))
#lapacksources = glob.glob(os.path.join('src/lapack', '*.c'))
#f2csources = glob.glob(os.path.join('src/f2c', '*.c'))

sources = splicingsources + lapacksources + f2csources

#new
include_dirs = [os.path.join(pysplicing_dir, 'include'),
                os.path.join(pysplicing_dir, 'src', 'lapack'),
                os.path.join(pysplicing_dir, 'src', 'f2c')]
#old
#include_dirs = ['include', 'src/lapack', 'src/f2c']

splicing_extension = Extension('pysplicing.pysplicing', sources, 
                               include_dirs=include_dirs,
                               define_macros=defines)


long_description = open('README').read()

if sys.version_info > (3, 0):
    options["use_2to3"] = True

setup(name = 'misopy',
      version = '0.1',
      description = 'Mixture of Isoforms model (MISO) for isoform quantitation using RNA-Seq',
      long_description = long_description,
#      license = 'MIT License',
      author = 'Yarden Katz,Gabor Csardi',
      author_email = 'gcsardi@stat.harvard.edu,yarden@mit.edu',
      url = 'http://genes.mit.edu/burgelab/miso/',
      ext_modules = [splicing_extension],
      # Tell distutils to look for pysplicing in the right directory
      package_dir = {'pysplicing': 'pysplicing/pysplicing'},
      packages = ['misopy', 'pysplicing'],
      # Required modules
      install_requires=[
#          "matplotlib >= 1.1.0",
          "matplotlib",
          "numpy >= 1.5.0",
          "scipy >= 0.9.0",
          "pysam >= 0.5.0"
          ],
      platforms = 'ALL',
      keywords = ['bioinformatics', 'sequence analysis',
                  'alternative splicing', 'RNA-Seq'],
      classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Programming Language :: C',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
        ]
      )

