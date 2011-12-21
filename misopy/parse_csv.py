##
## Helper functions for parsing csv files
##
## Yarden Katz, Sept 7, 2009
##

from numpy import *
from scipy import *
import time
import csv

def dictlist2csv(filename, dictlist, header_fields, delimiter='\t'):
    """
    Serialize a list of dictionaries into the output
    """
    str_header_fields = [str(f) for f in header_fields]
    header = "\t".join(str_header_fields) + '\n'
    output = open(filename, 'w')
    # write header to file
    output.write(header)
    for row in dictlist:
	row = "\t".join([str(row[field]) for field in header_fields]) + '\n'
	output.write(row)
    output.close()

def dictlist2dict(dictlist, header_name):
    """
    For the given dictlist, create a dictionary of each element keyed by
    the field in header_name.  Note that this assumes the header name is
    unique.
    """
    indexed_dict = {}
    for item in dictlist:
	indexed_dict[item[header_name]] = item
    return indexed_dict

def dictlist2array(dictlist, header_fields):
    """
    Convert a list of dictionaries into a numpy array, based on the given order of fields.
    """
    data_array = []
    for data_elt in dictlist:
	data_row = []
	for field in header_fields:
	    data_row.append(data_elt[field])
	data_array.append(data_row)
    return data_array

def csv2array(filename, skiprows=0, delimiter='\t', raw_header=False, missing=None, with_header=True):
    """
    Parse a file name into an array. Return the array and additional header lines. By default,
    parse the header lines into dictionaries, assuming the parameters are numeric,
    using 'parse_header'.
    """
    f = open(filename, 'r')
    skipped_rows = []
    for n in range(skiprows):
        header_line = f.readline().strip()
        if raw_header:
            skipped_rows.append(header_line)
        else:
            skipped_rows.append(parse_header(header_line))
    f.close()
    try:
        if missing:
            data = genfromtxt(filename, dtype=None, names=with_header,
                              deletechars='', skiprows=skiprows, missing=missing)
        else:
            if delimiter != '\t':
                data = genfromtxt(filename, dtype=None, names=with_header, delimiter=delimiter,
                                  deletechars='', skiprows=skiprows)
            else:
                data = genfromtxt(filename, dtype=None, names=with_header,
                                  deletechars='', skiprows=skiprows)
    except IOError as io_error:
        raise Exception, "IOError: %s. Filename: %s" %(io_error, filename)
    if data.ndim == 0:
	data = array([data.item()])
    return (data, skipped_rows)


def tryEval(s):
  try:
    return eval(s, {}, {})
  except:
    return s

def evalDict(d):
    for k, v in d.iteritems():
	d[k] = tryEval(v)
    return d


def get_header_fields(filename, delimiter='\t',
                      excel_tab=False):
    if excel_tab:
        f = open(filename, "rU")
    else:
        f = open(filename, "r")
    header_fields = f.readline().strip().split(delimiter)
    return header_fields


def file2dictlist(filename, delimiter='\t',
                  excel_tab=False):
    if excel_tab:
        f = open(filename, "rU")
        data = csv.DictReader(f, delimiter=delimiter,
                              quoting=csv.QUOTE_NONE,
                              dialect='excel')
        
    else:
        f = open(filename, "r")
        data = csv.DictReader(f, delimiter=delimiter,
                              quoting=csv.QUOTE_NONE)
    return data


def dictlist2file(dictrows, filename, fieldnames, delimiter='\t',
                  lineterminator='\n', extrasaction='ignore',
                  write_raw=False):
    out_f = open(filename, 'w')
    
    # Write out header
    if fieldnames != None:
        header = delimiter.join(fieldnames) + lineterminator
    else:
        header = dictrows[0].keys()
        header.sort()
    out_f.write(header)

    print "dictlist2file: serializing entries to %s" %(filename)

    t1 = time.time()
    if write_raw:
        for row in dictrows:
            out_f.write("%s%s" %(delimiter.join([row[name] for name in fieldnames]),
                                 lineterminator))
    else:
        # Write out dictionary
        data = csv.DictWriter(out_f, fieldnames,
                              delimiter=delimiter,
                              lineterminator=lineterminator,
                              extrasaction=extrasaction)
        for row in dictrows:
            data.writerow(row)
    out_f.close()
    t2 = time.time()
    print "dictlist2file: took %.2f seconds" %(t2 - t1)
    
    
def csv2dictlist_raw(filename, delimiter='\t'):
    f = open(filename)
    header_line = f.readline().strip()
    header_fields = header_line.split(delimiter)
    dictlist = []
    # convert data to list of dictionaries
    for line in f:
	values = map(tryEval, line.strip().split(delimiter))
	dictline = dict(zip(header_fields, values))
	dictlist.append(dictline)
    return (dictlist, header_fields)
    

def csv2dictlist(filename, raw_header=False, missing=None, delimiter=None,
                 with_header=True):
    """
    Parse a file into a list of dictionaries, where each dictionary has the values of that line
    keyed by the headers.
    """
    if not delimiter:
	data, header = csv2array(filename, raw_header=raw_header, missing=missing,
                                 with_header=with_header)
    else:
	data, header = csv2array(filename, raw_header=raw_header, missing=missing,
                                 delimiter=delimiter, with_header=with_header)
    header_line = open(filename).readline().strip()
    header_fields = header_line.split(delimiter)
    dictlist = []
    # convert data to list of dictionaries
    for values in data:
	dictline = dict(zip(header_fields, values))
	dictlist.append(dictline)
    return (dictlist, header_fields)

def find(val, values):
    """
    Find all instances of val in array. Return their indices.
    """
    indices = []
    values = list(values)
    n = 0
    for elt in values:
        if elt == val:
            indices.append(n)
        n += 1
    return indices

def parse_header(line, numeric_vals=True):
    """
    Parse a line of the form:
    
    #param=val\tparam=val\tparam=val...

    Return a dictionary of params: vals
    """
    line = line.strip()
    if line[0] == '#':
        line = line[1:]
    params = {}
    for pair in line.split('\t'):
        k, v = pair.split('=')
        if numeric_vals:
            params[k] = float(v)
        else:
            params[k] = v
    return params
