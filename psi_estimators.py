from scipy import *
from numpy import *

def curr_psi(NI, NE, skipped_exon_len, read_len, overhang_len, min_reads=0):
    """
    Compute the Psi value.                                        

    \Psi = DI / (DI + DE)                                              
    
    NI is: number of reads in upstream inclusion jxn + number of reads
    in downstream inclusion jxn + number of body reads in the skipped
    exon itself (for a skipped exon event, for example.)
  
    Note that NI + NE <= min_reads.  Return 'n/a' if this is not met.
    """
    if (NI + NE) < min_reads:
        raise Exception, "curr_psi: min_reads not satisfied."
    # First term of DI's denominator accounts for body reads, the
    # second for the two junctions (upstream, downstream)
    if skipped_exon_len - read_len < 0:
        raise Exception, "curr_psi: Exon too short for read length."
    DI = float(NI) / float((skipped_exon_len - read_len + 1) + 2*(read_len + 1 - (2*overhang_len)))
    # Denominator of DE accounts for the number of positions that our
    # read (minus the overhang) can be aligned to two adjacent
    # junctions
    DE = float(NE) / float(read_len + 1 - (2*overhang_len))
    if DI == 0 and DE == 0:
        raise Exception, "curr_psi undefined."
    psi = float(DI) / float(DI + DE)
    if psi > 1 or psi < 0:
        raise Exception, "Psi out of bounds [0, 1]."
    return psi

def get_psi_densities(NI, NE, skipped_exon_len, read_len, overhang_len):
    DI = float(NI) / float((skipped_exon_len - read_len + 1) + 2*(read_len + 1 - (2*overhang_len)))
    DE = float(NE) / float(read_len + 1 - (2*overhang_len))
    return DI, DE

def curr_psi_one_incjxn(NI, NE, skipped_exon_len, read_len, overhang_len, min_reads=0):
    """
    Compute the Psi value, using only one inclusion junction as the source of
    inclusion reads.

    \Psi = DI / (DI + DE)                                              
    
    NI is: number of reads in upstream inclusion jxn + number of reads
    in downstream inclusion jxn + number of body reads in the skipped
    exon itself (for a skipped exon event, for example.)
  
    Note that NI + NE <= min_reads.  Return 'n/a' if this is not met.
    """
    if (NI + NE) < min_reads:
       return 'n/a'

    DI = float(NI) / float(read_len + 1 - (2*overhang_len))
    # Denominator of DE accounts for the number of positions that our
    # read (minus the overhang) can be aligned to two adjacent
    # junctions
    DE = float(NE) / float(read_len + 1 - (2*overhang_len))
    ## Temporary
#    if DI == 0 and DE == 0:
#        return 0
    if DI == 0 and DE == 0:
	print NI, NE
        print "Estimator not defined for DI = 0, DE = 0."
#	return NaN
    psi = float(DI) / float(DI + DE)
    return psi

def curr_psi_bodyinc(NI, NE, skipped_exon_len, read_len, overhang_len, min_reads=0):
    """
    Compute the Psi value, using only body reads from skipped exon as source of
    inclusion reads.

    \Psi = DI / (DI + DE)                                              
    
    NI is: number of reads in upstream inclusion jxn + number of reads
    in downstream inclusion jxn + number of body reads in the skipped
    exon itself (for a skipped exon event, for example.)
  
    Note that NI + NE <= min_reads.  Return 'n/a' if this is not met.
    """
    if (NI + NE) < min_reads:
       return 'n/a'

    # First term of DI's denominator accounts for body reads, the
    # second for the two junctions (upstream, downstream)
    DI = float(NI) / float(skipped_exon_len - read_len + 1)
    # Denominator of DE accounts for the number of positions that our
    # read (minus the overhang) can be aligned to two adjacent
    # junctions
    DE = float(NE) / float(read_len + 1 - (2*overhang_len))
    if DI == 0 and DE == 0:
       return 0
    psi = float(DI) / float(DI + DE)
    return psi

def curr_psi_both_incjxns(NI, NE, skipped_exon_len, read_len, overhang_len, min_reads=0):
    """
    Compute the Psi value, using both inclusion junctions as source of inclusion reads.

    \Psi = DI / (DI + DE)                                              
    
    NI is: number of reads in upstream inclusion jxn + number of reads
    in downstream inclusion jxn + number of body reads in the skipped
    exon itself (for a skipped exon event, for example.)
  
    Note that NI + NE <= min_reads.  Return 'n/a' if this is not met.
    """
    if (NI + NE) < min_reads:
       return 'n/a'

    # First term of DI's denominator accounts for body reads, the
    # second for the two junctions (upstream, downstream)
    DI = float(NI) / float((read_len + 1 - (2*overhang_len)) + (read_len + 1 - (2*overhang_len)))
    # Denominator of DE accounts for the number of positions that our
    # read (minus the overhang) can be aligned to two adjacent
    # junctions
    DE = float(NE) / float(read_len + 1 - (2*overhang_len))
    if DI == 0 and DE == 0:
       return 0
    psi = float(DI) / float(DI + DE)
    return psi

def psi_bayes(n1, n2, nb, iso1_len, iso2_len, read_len, overhang_len, num_exons=3):
    """
    The MAP/MLE estimator for \Psi, based on a mixture model with uniform prior
    on \Psi.
    """
    p1 = 1/float((iso1_len - read_len + 1) - 4*(overhang_len-1))
    p2 = 1/float((iso2_len - read_len + 1) - 2*(overhang_len-1))
    if n1 == 0 and n2 == 0:
        raise Exception, "Psi Bayes undefined for n1, n2 == 0."
    psi = (n1*p1 + nb*p1 - 2*n1*p2 - n2*p2 - nb*p2 + sqrt(4*n1*p2*(n1*p1
          + n2*p1 + nb*p1 - n1*p2 - n2*p2 - nb*p2) + power((-n1*p1 - nb*p1 +
          2*n1*p2 + n2*p2 + nb*p2),2)))/(2*(n1*p1 + n2*p1 + nb*p1 - n1*p2 -
          n2*p2 - nb*p2))
    return psi

def transformed_psi_bayes(n1, n2, nb, iso1_len, iso2_len, read_len, overhang_len, num_exons=3):
    """
    Transformed MAP/MLE estimator for \Psi, taking into account lengths of
    isoforms.
    """
    if (n1 == 0 and n2 == 0) or (nb == 0):
        raise Exception, "transformed_psi_bayes not defined with: n1=%d, n2=%d, nb=%d" %(n1, n2, nb)
    psi_bayes_estimate = psi_bayes(n1, n2, nb, iso1_len, iso2_len, read_len, overhang_len, num_exons)
    ##
    ## Shouldn't take into account overhang, since it wasn't sampled with
    ## overhangs in mind
    ##
    num_positions_iso1 = (iso1_len - read_len + 1)
    num_positions_iso2 = (iso2_len - read_len + 1)
    trans_psi_bayes = (num_positions_iso2*psi_bayes_estimate)/(num_positions_iso1-num_positions_iso1*psi_bayes_estimate + num_positions_iso2*psi_bayes_estimate) 
    return trans_psi_bayes
