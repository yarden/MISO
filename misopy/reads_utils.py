##
## Utilities for manipulating reads and read assignments
##

from collections import defaultdict
from numpy import *

def count_aligned_reads(reads, paired_end=False):
    """
    Count the number of occurrences of each aligned read.
    Each read is a vector with K elements, where K is the
    number of isoforms, and each entry j is 1 or 0 depending
    on whether the read is consistent or inconsistent with
    the jth isoform.

    Returns the counts in a list of pairs, i.e.

    [(read_type1, counts), (read_type2, counts), ...]
    """
    counts_dict = defaultdict(int)

    for read in reads:
        if paired_end:
            hashable_read = tuple(map(int, read[0]))
        else:
            hashable_read = tuple(read)
        counts_dict[hashable_read] += 1

    # Sort results by keys for consistency
    keys = counts_dict.keys()
    keys.sort()

    counts = [(k, counts_dict[k]) for k in keys]
    
    return counts


def count_isoform_assignments(assignments):
    """
    Expects assignments to be an array of numbers.
    """
    num_isoforms = max(assignments)

    counts = [(iso_num, len(where(assignments == iso_num)[0])) \
              for iso_num in range(num_isoforms + 1)]

    return counts


def collapse_isoform_assignments(assignments):
    """
    Given a set of assignments like:

    [[0, 1],
     [1, 1]
     ...]

    Return a collapsed representation.
    """
    collapsed_assignments = defaultdict(int)
    for curr_assign in assignments:
        collapsed_assignments[curr_assign] += 1
    return collapsed_assignments
