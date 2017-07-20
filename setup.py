from setuptools import setup, Extension
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

# prefix directory for pysplicing module
pysplicing_dir = 'pysplicing'

splicingsources = glob.glob(os.path.join(pysplicing_dir, 'src', '*.c'))
lapacksources = glob.glob(os.path.join(pysplicing_dir, 'src', 'lapack', '*.c'))
f2csources = glob.glob(os.path.join(pysplicing_dir, 'src', 'f2c', '*.c'))

sources = splicingsources + lapacksources + f2csources

include_dirs = [os.path.join(pysplicing_dir, 'include'),
                os.path.join(pysplicing_dir, 'src', 'lapack'),
                os.path.join(pysplicing_dir, 'src', 'f2c')]

splicing_extension = Extension('pysplicing.pysplicing', sources, 
                               include_dirs=include_dirs,
                               define_macros=defines)

# Extract long description of MISO from README
long_description = open('README').read()

#if sys.version_info > (3, 0):
#    options["use_2to3"] = True

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
## Definition of the current version. This is defined here
## and then gets written to __init__.py in misopy.
##
MISO_VERSION = "0.5.4"

##
## Generate a __version__.py attribute
## for the module
##
with open("./misopy/__init__.py", "w") as version_out:
      version_out.write("__version__ = \"%s\"\n" %(MISO_VERSION))

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
      ext_modules = [splicing_extension],
      # Tell distutils to look for pysplicing in the right directory
      package_dir = {'pysplicing': 'pysplicing/pysplicing'},
      packages = ['misopy',
                  'misopy.sashimi_plot',
                  'misopy.sashimi_plot.plot_utils',
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
                    ('misopy/gff-events/mm9',
                     ['misopy/gff-events/mm9/SE.mm9.gff']),
                    ('misopy/gff-events/mm9/genes',
                     ['misopy/gff-events/mm9/genes/Atp2b1.mm9.gff']),
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

