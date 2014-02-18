##
## MISO into SQL database interface
##
import os
import sys
import time

import zipfile
import sqlite3
import shutil
import fnmatch
import glob

import misopy
import misopy.misc_utils as utils
import misopy.miso_utils as miso_utils


class MISODatabase:
    """
    Representation of a MISO SQLite database.
    """
    def __init__(self, db_fname):
        if not os.path.isfile(db_fname)
            raise Exception, "%s does not exist." %(db_fname)
        self.db_fname = db_fname
        

    def get_event_as_stream(self, chrom, event_name):
        """
        Get information for an event as a string
        stream.
        """
        pass


def miso_dir_to_db(miso_dir, output_fname):
    """
    Convert a MISO directory to an SQL db.
    """
    if not os.path.isdir(miso_dir):
        print "Error: %s not a directory." %(miso_dir)
        sys.exit(1)
    if os.path.isfile(output_fname):
        print "Error: %s already exists." %(output_fname)
        sys.exit(1)

        
def compress_miso_dir(self, dir_to_compress, output_filename):
    """
    Compress MISO directory into MySQL table using sqlite3.
    """
    print "Compressing %s into %s" %(dir_to_compress, output_filename)
    if not os.path.isdir(dir_to_compress):
        print "Error: %s not a directory, aborting." %(dir_to_compress)
        return None
    miso_filenames = glob.glob(os.path.join(dir_to_compress, "*.miso"))
    num_files = len(miso_filenames)
    print "  - %d files to compress" %(num_files)
    # Initialize the SQLite database
    if os.path.isfile(output_filename):
        print "Error: Compressed database %s already exists, aborting." \
              %(output_filename)
        return None
    conn = sqlite3.connect(output_filename)
    c = conn.cursor()
    # Create table for the current directory to compress
    table_name = os.path.basename(dir_to_compress)
    sql_create = \
        "CREATE TABLE %s " %(table_name) + \
        "(event_name text, psi_vals_and_scores text, header text)"
    c.execute(sql_create)
    for miso_fname in miso_filenames:
        miso_file_fields = load_miso_file_as_str(miso_fname)
        if miso_file_fields is None:
            print "Error: Cannot compress %s. Aborting." %(miso_fname)
            sys.exit(1)
        header, psi_vals_and_scores = miso_file_fields
        ######
        ###### TODO:
        ###### HANDLE COMPRESSED EVENT IDS HERE
        ######
        event_name = strip_miso_ext(os.path.basename(miso_fname))
        sql_insert = "INSERT INTO %s VALUES (?, ?, ?)" %(table_name)
        c.execute(sql_insert, (event_name,
                               psi_vals_and_scores,
                               header))
    # Commit changes and close the database
    conn.commit()
    conn.close()
    return output_filename
        

##
## Misc. helper functions
## 
def get_non_miso_files(filenames, miso_ext=".miso"):
    non_miso_files = []
    for fname in filenames:
        if os.path.basename(fname).endswith(miso_ext):
            non_miso_files.append(fname)
    return non_miso_files


def load_miso_file_as_str(miso_filename):
    """
    Load raw *.miso file as a set of strings to be inserted
    into an sqlite database.
    """
    if not os.path.isfile(miso_filename):
        print "Error: Cannot find %s" %(miso_filename)
        return None
    header = ""
    psi_vals_and_scores = ""
    with open(miso_filename) as miso_file:
        # Read the header, consisting of two lines
        for n in range(2):
            header += miso_file.readline()
        for line in miso_file:
            psi_vals_and_scores += line
    return header, psi_vals_and_scores
                

def strip_miso_ext(filename):
    if filename.endswith(".miso"):
        return filename[0:-5]
    return filename
