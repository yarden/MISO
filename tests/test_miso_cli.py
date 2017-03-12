
import os
import sys
import shutil

from nose.tools import *

def miso_cli(sams = None, gff = None):

    import misopy
    miso_dir = misopy.__path__[0]

    # Data files to test directory, if the output directory exists,
    # then we assume that it has these files already.
    if gff is None:
        test_gff = os.path.join(miso_dir, "..", "tests", "test.gff")
    if sams is None:
        sams = [ os.path.join(miso_dir, "test-data", "sam-data",
                             "c2c12.Atp2b1.sam") ]

    dirpath = 'output'
    curdir = os.getcwd()
    try:
        os.mkdir(dirpath)
    except OSError:
        None

    try:
        os.chdir(dirpath)
        shutil.copyfile(test_gff, os.path.basename(test_gff))
        [ shutil.copyfile(sam, os.path.basename(sam)) for sam in sams ]

        # SAM to BAM
        from misopy import sam_to_bam
        def tobam(sam):
            sam_to_bam.main(["--convert", os.path.basename(sam), "my_sample"])
            return os.path.join(
                "my_sample",
                os.path.splitext(os.path.basename(sam))[0] + ".sorted.bam")

        bams = [ tobam(sam) for sam in sams ]

        # Index GFF
        from misopy import index_gff
        index_gff.main(["--index", os.path.basename(test_gff), "indexed"])

        # Run MISO
        from misopy import miso
        miso_args = ["--run", "indexed"] + bams + \
                    ["--output-dir", "output", "--read-len", "36"]
        miso.main(miso_args)

        # Summarize
        from misopy import summarize_miso
        summarize_miso.main(["--summarize-samples", "output", "."])
    finally:
        os.chdir(curdir)

def test_miso_cli():
    miso_cli()

@raises(SystemExit)
def test_error_if_multiple_bam_files():
    import misopy
    miso_dir = misopy.__path__[0]
    test_sam = os.path.join(miso_dir, "test-data", "sam-data",
                            "c2c12.Atp2b1.sam")
    test_sam2 = os.path.join(os.path.dirname(test_sam),
                             "c2c12.Atp2b1-2.sam")
    try:
        shutil.copyfile(test_sam, test_sam2)
        miso_cli(sams = [test_sam, test_sam2])
    finally:
        if os.path.isfile:
            os.remove(test_sam2)
