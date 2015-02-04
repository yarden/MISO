
import os
import sys

import misopy

misopy_dir = os.path.dirname(misopy.__file__)
miso_dir = os.path.dirname(misopy_dir)
pysplicing_dir = os.path.join(miso_dir, "pysplicing")

sys.path.insert(0, miso_dir)
sys.path.insert(0, misopy_dir)
sys.path.insert(0, pysplicing_dir)

os.environ["PYTHONPATH"] = os.pathsep.join(sys.path)

import pysplicing
