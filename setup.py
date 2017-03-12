from setuptools import setup, Extension
import distutils.ccompiler
import glob
import os
import sys
import shutil

if not os.path.exists("temp"):
    os.makedirs("temp")

## Test for functions, with a hack to suppress compiler warnings.
os.chdir('temp')
cc = distutils.ccompiler.new_compiler()
cc.define_macro('main', 'int (main)')
defines = []
if cc.has_function('float a = rintf(1.0); return rand',
                   includes=['math.h', 'stdlib.h'], libraries=['m']):
    defines.append(('HAVE_RINTF', '1'))

if cc.has_function('int a = isfinite(1.0); return rand',
                   includes=['math.h', 'stdlib.h'], libraries=['m']):
    defines.append(('HAVE_ISFINITE', '1'))

if cc.has_function('double a = expm1(1.0); return rand',
                   includes=['math.h', 'stdlib.h'], libraries=['m']):
    defines.append(('HAVE_EXPM1', '1'))

if cc.has_function('double a = rint(1.0); return rand',
                   includes=['math.h', 'stdlib.h'], libraries=['m']):
    defines.append(('HAVE_RINT', '1'))

if cc.has_function('double log2(double); double a = log2(1.0); return rand',
                   includes=['math.h', 'stdlib.h'], libraries=['m']):
    defines.append(('HAVE_LOG2', '1'))

if cc.has_function('long double a = logbl(1.0); return rand',
                   includes=['math.h', 'stdlib.h'], libraries=['m']):
    defines.append(('HAVE_LOGBL', '1'))

if cc.has_function('snprintf(0, 0, ""); return rand',
                   includes=['stdio.h', 'stdlib.h']):
    defines.append(('HAVE_SNPRINTF', '1'))

if cc.has_function('double a = log1p(1.0); return rand',
                   includes=['math.h', 'stdlib.h'], libraries=['m']):
    defines.append(('HAVE_LOG1P', '1'))

if cc.has_function('double round(double) ; double a = round(1.0); return rand', 
                   includes=['math.h', 'stdlib.h'], libraries=['m']):
    defines.append(('HAVE_ROUND', '1'))

if cc.has_function('double fmin(double, double); double a = fmin(1.0,0.0); return rand', 
                   includes=['math.h', 'stdlib.h'], libraries=['m']):
    defines.append(('HAVE_FMIN', '1'))

os.chdir('..')
shutil.rmtree("temp")

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
MISO_VERSION = "0.5.3"

## This is our implementation of dependency tracking,
## to avoid rebuilding the C files continuesly.
## Idea is from http://stackoverflow.com/questions/11013851

def my_compile(self, sources, output_dir=None, macros=None,
               include_dirs=None, debug=0, extra_preargs=None,
               extra_postargs=None, depends=None):

    macros, objects, extra_postargs, pp_opts, build = \
        self._setup_compile(output_dir, macros, include_dirs, sources,
                            depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)

    # -------------------------------------------------------------
    # This is where we intervine. (The rest is the original
    # CCompiler.compile.) We check which file in objects really
    # needs a rebuild, and only keep the the ones that do.

    # We need to create a temporary Makefile, with all the
    # dependency data, collected from a previous compilation
    import tempfile
    import re
    import subprocess

    makefile_os = tempfile.mkstemp()
    makefile = os.fdopen(makefile_os[0], "w")

    makefile.write('objects: ' + ' '.join(objects) + '\n' +
                   '.PHONY: objects\n\n')
    depfiles = ' '.join([ re.sub("\.o$", ".d", of) for of in objects ])
    makefile.write('-include ' + depfiles + '\n\n')

    for obj in objects:
        try:
            src, ext = build[obj]
        except KeyError:
            continue
        makefile.write(obj + ': ' + src + '\n' +
                       '\t@echo $@\n\n')

    makefile.close()

    ood_objects = subprocess.check_output(['make', 'objects', '-f',
                                           makefile_os[1]]).strip().split('\n')

    # We need to add compiler flags to create the dependency files
    cc_args += ['-MMD', '-MP']

    # -------------------------------------------------------------

    for obj in ood_objects:
        try:
            src, ext = build[obj]
        except KeyError:
            continue
        self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)

        # Return *all* object filenames, not just the ones we just built.
    return objects

try:
    if os.environ['GABOR'] == 'yes':
        distutils.ccompiler.CCompiler.compile = my_compile
except:
    None

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
      author_email = 'yarden@mit.edu,csardi.gabor@gmail.com',
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

