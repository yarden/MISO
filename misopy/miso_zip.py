##
## Compress/uncompress directories containing MISO output.
##
## Reduces the number of raw output files (*.miso files) by
## storing them in a database, and then compressing using
## zip compression to reduce the space footprint.
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


class MISOCompressor:
    """
    Compressor/uncompressor of MISO output-containing directories.
    """
    def __init__(self):
        self.input_dir = None
        self.output_dir = None
        # Extension for compressed directory that contains
        # any sort of MISO output within it
        self.comp_ext = ".misozip"
        # Extension for compressed chromosome directories
        self.miso_db_ext = ".miso_db"
        

    def compress(self, output_filename, miso_dirnames):
        """
        Takes a MISO input directory and compresses it
        into 'output_filename'.
        """
        if os.path.isfile(output_filename):
            print "Error: %s already exists. Please delete to overwrite." \
                  %(output_filename)
        output_dir = "%s%s" %(output_filename, self.miso_db_ext)
        if os.path.isdir(output_dir):
            print "Error: Intermediate compressed directory %s " \
                  "exists. Please delete to overwrite." %(output_dir)
            sys.exit(1)
        for miso_dirname in miso_dirnames:
            print "Processing: %s" %(miso_dirname)
            if not os.path.isdir(miso_dirname):
                print "Error: %s not a directory." %(miso_dirname)
                sys.exit(1)
            if os.path.isfile(output_filename):
                print "Output file %s already exists, aborting. " \
                      "Please delete the file if you want " \
                      "compression to run."
                sys.exit(1)
            self.miso_dirs_to_compress = []
            print "Copying source directory tree.."
            shutil.copytree(miso_dirname, output_dir,
                            ignore=self.collect_miso_dirs)
            for dir_to_compress in self.miso_dirs_to_compress:
                rel_path = os.path.relpath(dir_to_compress, miso_dirname)
                comp_path = os.path.join(output_dir, rel_path)
                # Remove the place holder directory
                os.rmdir(comp_path)
                # Create the zip file
                comp_path = "%s%s" %(comp_path, self.miso_db_ext)
                self.compress_miso_dir(dir_to_compress, comp_path)
        # Zip directory using conventional zip
        print "Zipping compressed directory with standard zip..."
        t1 = time.time()
        zipper(output_dir, output_filename)
        print "Deleting intermediate directory: %s" %(output_dir)
        shutil.rmtree(output_dir)
        t2 = time.time()
        print "  - Standard zipping took %.2f minutes." \
              %((t2 - t1)/60.)
        print "To access the SQLite representation of raw MISO output "
        print "(*.miso) files, simply unzip with the .miso_zip file "
        print "with standard unzip utility:\n"
        print "  unzip %s" %(output_filename)
        

    def uncompress(self, compressed_filename, output_dir):
        """
        Takes a compressed MISO file 'compressed_filename' and
        uncompresses it into 'output_dir'.
        """
        if not os.path.isfile(compressed_filename):
            print "Error: Cannot find %s, aborting." \
                  %(compressed_filename)
        if not os.path.basename(compressed_filename).endswith(self.comp_ext):
            print "Warning: %s does not end in %s. Are you sure it is " \
                  "a file compressed by miso_zip.py?" \
                  %(compressed_filename, self.comp_ext)
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        print "Uncompressing %s into %s" %(compressed_filename,
                                           output_dir)
        # First unzip the file using conventional unzip
        unzipped_files = unzipper(compressed_filename, output_dir)
        self.uncompressed = {}
        for unzipped_fname in unzipped_files:
            self.uncompress_miso_dir(unzipped_fname, output_dir)
        # Remove the original .zip
        if os.path.isfile(compressed_filename):
            print "Removing the compressed file %s" %(compressed_filename)
            if os.path.isfile(compressed_filename):
                os.remove(compressed_filename)
            

    def uncompress_miso_dir(self, compressed_filename, output_dir):
        """
        Uncompress a MISO compressed dir (*.miso_db) file. Replaces
        file 'name.miso_db' file with a directory named 'name' that
        contains a .miso file for each event contained in 'name.miso_db'
        """
        if os.path.isdir(compressed_filename):
            for root, dirs, files in os.walk(compressed_filename):
                for filename in files:
                    if self.is_miso_db_file(filename):
                        fname_to_uncompress = os.path.join(root, filename)
                        self.uncompress_miso_file(fname_to_compress,
                                                  output_dir)
        elif self.is_miso_db_file(compressed_filename):
            self.uncompress_miso_file(compressed_filename,
                                      output_dir)


    def uncompress_miso_file(self, compressed_filename, output_dir):
        print "Uncompressing %s.." %(compressed_filename)
        compressed_filename = os.path.join(output_dir,
                                           compressed_filename)
        table_name = self.get_table_name_from_file(compressed_filename)
        if table_name is None:
            print "Error: Cannot retrieve name of MISO db file %s" \
                  %(compressed_filename)
            return None
        conn = sqlite3.connect(compressed_filename)
        c = conn.cursor()
        results = c.execute("SELECT * from %s" %(table_name))
        # Create an uncompressed directory that contains the
        # contents of the current *.miso_db file
        uncompressed_dir = \
            os.path.join(os.path.dirname(compressed_filename),
                         table_name)
        print "Making uncompressed directory %s" %(uncompressed_dir)
        if not os.path.isdir(uncompressed_dir):
            os.makedirs(uncompressed_dir)
        for row in results:
            event_name, psi_vals_and_scores, header = row
            # Create *.miso file for current event
            miso_filename = os.path.join(uncompressed_dir,
                                         "%s.miso" %(event_name))
            with open(miso_filename, "w") as miso_file:
                # Write header
                miso_file.write(header)
                # Write psi vals and their log scores
                miso_file.write(psi_vals_and_scores)
        # Remove the MISO db file once done uncompressing it
        print "Deleting: %s" %(compressed_filename)
        os.remove(compressed_filename)
        conn.close()


    def get_table_name_from_file(self, db_filename):
        """
        Return the name of a file.
        """
        base_name = os.path.basename(db_filename)
        if base_name.endswith(self.miso_db_ext):
            # Strip extension
            return base_name[0:-1*len(self.miso_db_ext)]
        return None
        

    def is_miso_db_file(self, compressed_filename):
        if os.path.isdir(compressed_filename):
            return False
        if os.path.basename(compressed_filename).endswith(self.miso_db_ext):
            return True
        return False
            

    def collect_miso_dirs(self, path, dirnames):
        """
        Collect raw MISO output directories
        """
        if fnmatch.filter(dirnames, "*.miso"):
            self.miso_dirs_to_compress.append(path)
            return dirnames
        return []


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
                

def compress_miso(output_filename, input_dirs,
                  comp_ext=".misozip"):
    """
    Compress a directory containing MISO files.

    Traverse directories, one by one, and look for directories
    that contain 
    """
    output_filename = utils.pathify(output_filename)
    for input_dir in input_dirs:
        if not os.path.isdir(input_dir):
            print "Error: Cannot find directory %s" %(input_dir)
            sys.exit(1)
    if not os.path.basename(output_filename).endswith(comp_ext):
        print "Error: Compressed output filename must end in %s" \
              %(comp_ext)
        sys.exit(1)
    if os.path.isfile(output_filename):
        print "Error: Output filename exists. Please delete %s to overwrite." \
              %(output_filename)
        sys.exit(1)
    t1 = time.time()
    miso_comp = MISOCompressor()
    miso_comp.compress(output_filename, input_dirs)
    t2 = time.time()
    print "Compression took %.2f minutes." %((t2 - t1)/60.)
    

def uncompress_miso(compressed_filename, output_dir):
    """
    Uncompress MISO directory.
    """
    if not os.path.isfile(compressed_filename):
        print "Error: Zip file %s is not found." \
              %(compressed_filename)
        sys.exit(1)
    t1 = time.time()
    miso_comp = MISOCompressor()
    miso_comp.uncompress(compressed_filename, output_dir)
    t2 = time.time()
    print "Uncompression took %.2f minutes." %((t2 - t1)/60.)


def strip_miso_ext(filename):
    if filename.endswith(".miso"):
        return filename[0:-5]
    return filename


def zipper(dir, zip_file):
    """
    Zip a directory 'dir' recursively, saving result in
    'zip_file'.
    
    by Corey Goldberg.
    """
    zip = zipfile.ZipFile(zip_file, 'w',
                          compression=zipfile.ZIP_DEFLATED)
    root_len = len(os.path.abspath(dir))
    for root, dirs, files in os.walk(dir):
        archive_root = os.path.abspath(root)[root_len:]
        for f in files:
            fullpath = os.path.join(root, f)
            archive_name = os.path.join(archive_root, f)
            zip.write(fullpath, archive_name, zipfile.ZIP_DEFLATED)
    zip.close()
    return zip_file


def unzipper(zip_file, outdir):
    """
    Unzip a given 'zip_file' into the output directory 'outdir'.
    
    Return the names of files in the archive.
    """
    zf = zipfile.ZipFile(zip_file, "r")
    filenames = zf.namelist()
    zf.extractall(outdir)
    return filenames


def greeting():
    print "Compress/uncompress MISO output. Usage:\n"
    print "To compress a directory containing MISO files \'inputdir\', use: "
    print "   miso_zip.py --compress outputfile.misozip inputdir"
    print "To uncompress back into a directory \'outputdir\', use: "
    print "   miso_zip --uncompress outputfile.misozip outputdir"
    print "\nNote: compressed filename must end in \'.misozip\'"
    

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--compress", dest="compress", nargs=2, default=None,
                      help="Compress a directory containing MISO output. "
                      "Takes as arguments: (1) the output filename of the "
                      "compressed file, (2) a comma-separated list of "
                      "directory names to be compressed. "
                      "Example: --compress output.misozip dirname1,dirname2")
    parser.add_option("--uncompress", dest="uncompress", nargs=2, default=None,
                      help="Uncompress a file generated by compress_miso. "
                      "Takes as arguments: (1) the filename to be "
                      "uncompressed, and (2) the directory to place the "
                      "uncompressed representation into. "
                      "Example: --uncompress output.misozip outputdir")
    (options, args) = parser.parse_args()

    if (options.compress is None) and (options.uncompress is None):
        greeting()
        sys.exit(1)
    elif (options.compress is not None) and (options.uncompress is not None):
        # Can't be given both.
        greeting()
        print "Error: Cannot process --compress and --uncompress at same time."
        sys.exit(1)

    if options.compress is not None:
        output_filename = utils.pathify(options.compress[0])
        input_dirs = [utils.pathify(d) \
                      for d in options.compress[1].split(",")]
        compress_miso(output_filename, input_dirs)

    if options.uncompress is not None:
        compressed_filename = utils.pathify(options.uncompress[0])
        output_dir = utils.pathify(options.uncompress[1])
        uncompress_miso(compressed_filename, output_dir)


if __name__ == "__main__":
    main()
    
