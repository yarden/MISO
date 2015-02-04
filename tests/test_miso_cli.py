
import os
import sys
import tempfile
import shutil

from nose.tools import *

import setup_miso

def miso_cli(sams = None, gff = None):

    miso_dir = setup_miso.miso_dir

    # Data files in a temporary directory
    if gff is None:
        test_gff = os.path.join(miso_dir, "misopy", "gff-events", "mm9",
                                "SE.mm9.gff")
    if sams is None:
        sams = [ os.path.join(miso_dir, "misopy", "test-data", "sam-data",
                             "c2c12.Atp2b1.sam") ]

    dirpath = tempfile.mkdtemp()
    os.chdir(dirpath)
    shutil.copyfile(test_gff, os.path.basename(test_gff))
    [ shutil.copyfile(sam, os.path.basename(sam)) for sam in sams ]

    # SAM to BAM
    import sam_to_bam
    def tobam(sam):
        sam_to_bam.main(["--convert", os.path.basename(sam), "my_sample"])
        return os.path.join(
            "my_sample",
            os.path.splitext(os.path.basename(sam))[0] + ".sorted.bam")

    bams = [ tobam(sam) for sam in sams ]

    # Index GFF
    import index_gff
    index_gff.main(["--index", os.path.basename(test_gff), "indexed"])

    # Run MISO
    import miso
    miso_args = ["--run", "indexed"] + bams + \
                ["--output-dir", "output", "--read-len", "36"]
    miso.main(miso_args)

    # Summarize
    import summarize_miso
    summarize_miso.main(["--summarize-samples", "output", "."])

    # Clean up
    shutil.rmtree(dirpath)

def test_miso_cli():
    miso_cli()

@raises(SystemExit)
def test_error_if_multiple_bam_files():
    miso_dir = setup_miso.miso_dir
    test_sam = os.path.join(miso_dir, "misopy", "test-data", "sam-data",
                            "c2c12.Atp2b1.sam")
    test_sam2 = os.path.dirname(test_sam) + "c2c12.Atp2b2-2.sam"
    shutil.copyfile(test_sam, test_sam2)

    miso_cli(sams = [test_sam, test_sam2])
