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
import StringIO

import misopy
import misopy.misc_utils as misc_utils
import misopy.miso_utils as miso_utils

# File extension for MISO SQLite databases
MISO_DB_EXT = ".miso_db"


class MISODatabase:
    """
    Representation of a MISO SQLite database.
    """
    def __init__(self, db_fname, comp_to_uncomp=None):
        self.comp_to_uncomp = comp_to_uncomp
        ####
        #### TODO MAKE MAPPING FROM UNCOMPRESSED TO COMPRESSED
        #### GENES HERE
        ####
        self.uncomp_to_comp = None
        if self.comp_to_uncomp is not None:
            # Make mapping from uncompressed to compressed IDs
            self.uncomp_to_comp = misc_utils.inv_dict(self.comp_to_uncomp)
        if not os.path.isfile(db_fname):
            raise Exception, "%s does not exist." %(db_fname)
        self.db_fname = db_fname
        self.table_name = get_table_name_from_file(self.db_fname)
        if self.table_name is None:
            print "Error: Cannot retrieve name of MISO db file %s" \
                  %(self.db_fname)
            return None
        self.conn = sqlite3.connect(self.db_fname)


    def get_event_data_as_stream(self, event_name):
        """
        Get data for given event. If there's no data, return None.
        """
        # The name of the event as stored in database (if using
        # compressed event IDs, this would be a misocomp ID)
        event_to_query = event_name
        # If we have a compressed event mapping, use it
        if self.uncomp_to_comp is not None:
            if misc_utils.is_compressed_name(event_name):
                # If event name is already compressed, use it
                event_to_query = event_name
                # Event to return is uncompressed name
                if event_name not in self.comp_to_uncomp:
                    print "WARNING: Failed to find uncompressed ID for %s" \
                          %(event_name)
                    return None
            else:
                # Event name is not compressed
                if event_name not in self.uncomp_to_comp:
                    print "WARNING: Failed to find compressed ID for %s" \
                          %(event_name)
                    return None
                event_to_query = self.uncomp_to_comp[event_name]
        c = self.conn.cursor()
        results = \
          c.execute("SELECT * from %s WHERE event_name=\'%s\'" \
                    %(self.table_name,
                      event_to_query))
        rows = results.fetchall()
        if len(rows) == 0:
            # Event not found
            return None
        if len(rows) > 1:
            raise Exception, \
              "More than one entry for event %s" %(event_to_query)
        event_name, psi_vals_and_scores, header = rows[0]
        # If we're given a mapping to compressed events,
        # return the event name as the *uncompressed* event
        # name
        event_data = "%s\n%s\n" %(header,
                                  psi_vals_and_scores)
        print "EVENT DATA: ", event_data
        event_stream = StringIO.StringIO(event_data)
        return event_stream


    def get_event_data_as_string(self, event_name):
        data = self.get_event_data_as_stream(event_name).read()
        return data


    def get_all_events(self):
        """
        Get all events from table. Return iterator of results.
        """
        c = self.conn.cursor()
        results = c.execute("SELECT * from %s" %(self.table_name))
        return results

    
    def get_all_event_names(self):
        """
        Return all event names
        """
        event_names = []
        for result in self.get_all_events():
            event_names.append(result[0])
        return event_names

        
def miso_dir_to_db(dir_to_compress, output_filename):
    """
    Convert MISO directory into MySQL table using sqlite3.
    """
    print "Converting MISO directory into database" 
    print "  - MISO dir: %s" %(dir_to_compress)
    print "  - Output file: %s" %(output_filename)
    if not os.path.isdir(dir_to_compress):
        print "Error: %s not a directory, aborting." %(dir_to_compress)
        sys.exit(1)
    miso_filenames = glob.glob(os.path.join(dir_to_compress, "*.miso"))
    num_files = len(miso_filenames)
    print "  - %d files to compress" %(num_files)
    # Initialize the SQLite database
    if os.path.isfile(output_filename):
        print "Error: Database %s already exists, aborting." \
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
            return None
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


def get_table_name_from_file(db_filename):
    """
    Return the name of a file.
    """
    base_name = os.path.basename(db_filename)
    if base_name.endswith(MISO_DB_EXT):
        # Strip extension
        return base_name[0:-1*len(MISO_DB_EXT)]
    return None


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


def is_miso_db_fname(fname):
    """
    Return True if it's a MISO db filename, like
    chrX.miso_db
    """
    return fname.endswith(MISO_DB_EXT)


def is_miso_unpacked_dir(dirname):
    """
    Return True if it's a MISO directory that is unpacked,
    i.e. a directory containing *.miso text files immediately
    within it (so having these files in a subdir of that dir
    will not count.)
    """
    if not os.path.isdir(dirname):
        return False
    matches = fnmatch.filter(os.listdir(dirname), "*.miso")
    if len(matches) != 0:
        return True
    return False
                

def strip_miso_ext(filename):
    if filename.endswith(".miso"):
        return filename[0:-5]
    return filename
