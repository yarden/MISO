
[![Build Status](https://travis-ci.org/gaborcsardi/MISO.svg?branch=fastmiso)](https://travis-ci.org/gaborcsardi/MISO)

MISO (Mixture of Isoforms) is a probabilistic framework that quantitates the expression level of alternatively spliced genes from RNA-Seq data, and identifies differentially regulated isoforms or exons across samples. By modeling the generative process by which reads are produced from isoforms in RNA-Seq, the MISO model uses Bayesian inference to compute the probability that a read originated from a particular isoform.

MISO uses the inferred assignment of reads to isoforms to quantitate the abundances of the underlying set of alternative mRNA isoforms. Confidence intervals over estimates can be obtained, which quantify the reliability of the estimates. 

Please see http://genes.mit.edu/burgelab/miso/docs/ for MISO manual and documentation.
