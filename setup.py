from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import glob
import os
import sys

# Extract long description of MISO from README
long_description = open('README').read()

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
MISO_VERSION = "0.4.8python"

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
      packages = ['misopy',
                  'misopy.sashimi_plot',
                  'misopy.sashimi_plot.plot_utils'],
      # Cython extensions
      cmdclass = {'build_ext': build_ext},
      ext_modules = [Extension("misopy.miso_scores", ["misopy/miso_scores.pyx"])],
      # distutils always uses forward slashes
      scripts = ['misopy/module_availability.py',
                 'misopy/index_gff.py',
                 'misopy/sam_to_bam.py',
                 'misopy/run_events_analysis.py',
                 'misopy/run_miso.py',
                 'misopy/exon_utils.py',
                 'misopy/pe_utils.py',
                 'misopy/filter_events.py',
                 'misopy/miso_zip.py',
                 # sashimi_plot scripts
                 'misopy/sashimi_plot/plot.py'],
      data_files = [('misopy/settings',
                     ['misopy/settings/miso_settings.txt',
                      'misopy/sashimi_plot/settings/sashimi_plot_settings.txt']),
                    ('misopy/test-data/sam-data',
                     ['misopy/test-data/sam-data/c2c12.Atp2b1.sam']),
                    ('misopy/gff-events/mm9',
                     ['misopy/gff-events/mm9/SE.mm9.gff']),
                    ('misopy/gff-events/mm9/genes',
                     ['misopy/gff-events/mm9/genes/Atp2b1.mm9.gff']),
                    ('misopy/gff-events/mm9/events',
                     ['misopy/gff-events/mm9/events/test_event.gff']),                     
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
          "cython",
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

