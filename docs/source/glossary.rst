.. include:: <isogrk3.txt>

.. contents::


Glossary of terms related to MISO
=================================

Terms used in the `MISO manual <index.html>`_. 

Annotations and GFF files
-------------------------

* **Exon-centric**: An annotation of alternative events in the genome that is based on inclusion/exclusion of a particular exon in a transcript. For example, an exon-centric annotation of an alternatively skipped exon would contain two isoforms, one containing the skipped exon and its two flanking exon, and another isoform containing only the two flanking exons. This "exon-centric" annotation does not incorporate other exons in the gene, and so |Psi| values obtained from this annotation correspond only to the inclusion of the alternative exon relative to its two annotated flanking exons, without considering any other parts of the gene's isoforms.

* **Isoform-centric**: Unlike exon-centric annotations, in isoform-centric annotations each whole isoform of a gene is annotated and used as input to MISO. |Psi| values obtained this way are vectors, each entry corresponding to the percent inclusion of a whole isoform in the annotated gene.

Read alignments and BAM files
-----------------------------

* **Paired-end versus single-end**: In paired-end sequencing, both ends of a cluster on a flow cell are sequenced. Each mate is guaranteed to have originated from the same molecule. In single-end sequencing, only one end of a molecule is sequenced. MISO supports both paired-end and single-end data. All paired-end data can be run in MISO as single-end by simply omitting the ``--paired-end`` parameter. In that case, MISO will treat each mate of a pair independently.

* **Properly paired read pair**: This term applies only to paired-end data, and refers to read pairs where both mates are mapped in a way that makes sense given the strandedness of the RNA-Seq protocol and the alignments of the individual mates. When MISO maps read pairs to event annotations in paired-end mode, it only considers properly paired reads. If the mates maps to distinct chromosomes, then the read pair will not be considered properly paired. Similarly, if one mate maps in opposite orientation to what is expected given the strandedness setting, it will not be considered properly paired. Finally, if one mate maps to within the boundaries of an annotated event but the other does not, the read pair will not be considered (though if such cases are common, one can use MISO in single-end mode.) MISO will generally look for the BAM flag that encodes whether a read pair is properly paired or not. Otherwise, it pairs mates together from a BAM file using their read IDs. 

* **Overhang**: Overhang applied to splice junction reads, refers to the minimum number of bases covered by the read on any of the exons involved in the junction. For example, if a junction read of length 30 is aligned to the border of two exons with 10 bases covered on one exon and 20 bases covered on the other exon, the overhang is defined to be 10 (the smallest of the two numbers.) For single-end reads, requiring a considerable overhang like 4 or more helps filter alignments that appear as junction reads but are simply artifacts of sequencing errors and/or alignment errors. Overhang is not defined for paired-end reads.


Inference terms
---------------

There are a number of technical parameter settings related to Markov chain Monte Carlo inference (MCMC),
which the MISO engine is based on. In virtually all cases, users never have to mind or alter these
settings, but they are explained here for completeness. These are configurable from the MISO settings file.

* **Number of (MCMC) chains**: The number of independent MCMC chains used by MISO when performing inferences.
The default number is 6 which is considered a conservative setting for the problem. High chain numbers like 6
prevent MISO from getting stuck in suboptimal |Psi| values.

* **Lag**: The number of MCMC samples to skip over when computing the posterior distribution over |Psi|. The default is 10. High settings of this parameter can, in some cases, prevent autocorrelations between MCMC samples.

* **Burn-in**: The number of initial MCMC samples to exclude when computing the posterior distribution over |Psi|. Large settings of burn-in can prevent generation of posterior distributions over |Psi| that are closely correlated to the initial random setting of |Psi| used by the sampler.


.. _MISO: index.html
.. _A primer on probabilistic inference: http://cocosci.berkeley.edu/tom/papers/tutorial.pdf

.. _distribute: http://packages.python.org/distribute/
.. _pip: http://pypi.python.org/pypi/pip
.. _Enthought Python Distribution: http://www.enthought.com/products/epd.php
.. _EPD: http://www.enthought.com/products/epd.php
.. _Superpack: http://fonnesbeck.github.com/ScipySuperpack/
.. _UCSC Genome Browser: http://genome.ucsc.edu/
.. _SQLite database: http://www.sqlite.org/
.. _BioMart: http://www.ensembl.org/biomart/martview 
.. _A Practical Course in Bayesian Graphical Modeling: http://www.socsci.uci.edu/~mdlee/bgm.html
.. _Python 2.6: http://www.python.org
.. _Numpy: http://www.numpy.org
.. _Scipy: http://www.scipy.org
.. _matplotlib: http://matplotlib.sourceforge.net/
.. _simplejson: http://code.google.com/p/simplejson
.. _pygsl: http://pygsl.sourceforge.net
.. _GSL: http://www.gnu.org/software/gsl
.. _samtools: http://samtools.sourceforge.net/
.. _pysam: http://code.google.com/p/pysam/
.. _Spliced Alignment/MAP (SAM): http://samtools.sourceforge.net/SAM1.pdf
.. _SAM: http://samtools.sourceforge.net/SAM1.pdf
.. _GFF: http://www.sequenceontology.org/gff3.shtml
.. _GTF2GFF: http://www.sequenceontology.org/resources/converter_readme.html
.. _gtf2gff3.pl: http://genes.mit.edu/burgelab/miso/scripts/gtf2gff3.pl
.. _ucsc_table2gff3.pl: http://genes.mit.edu/burgelab/miso/scripts/ucsc_table2gff3.pl
.. _genePredToGtf: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
.. _biotoolbox: http://code.google.com/p/biotoolbox/
.. _Drosophila melanogaster alternative events (modENCODE): http://genes.mit.edu/burgelab/miso/annotations/modENCODE_alt_events.zip
.. _Events to genes mapping for mm9: http://genes.mit.edu/burgelab/miso/annotations/mm9_events_to_genes.zip
.. _Events to genes mapping for hg18: http://genes.mit.edu/burgelab/miso/annotations/hg18_events_to_genes.zip
.. _Events to genes mapping for hg19: http://genes.mit.edu/burgelab/miso/annotations/hg19_events_to_genes.zip
.. _mm9 ensGene GFF annotation: http://genes.mit.edu/burgelab/miso/annotations/ucsc_tables/mm9/ensGene.gff3
.. _mm10 ensGene GFF annotation: http://genes.mit.edu/burgelab/miso/annotations/ucsc_tables/mm10/ensGene.gff3
.. _hg18 ensGene GFF annotation: http://genes.mit.edu/burgelab/miso/annotations/ucsc_tables/hg18/ensGene.gff3
.. _hg19 ensGene GFF annotation: http://genes.mit.edu/burgelab/miso/annotations/ucsc_tables/hg19/ensGene.gff3
.. _Mouse genome (mm9) alternative events: http://genes.mit.edu/burgelab/miso/annotations/mm9_alt_events.zip
.. _Human genome (hg18) alternative events: http://genes.mit.edu/burgelab/miso/annotations/hg18_alt_events.zip
.. _Human genome (hg19) alternative events: http://genes.mit.edu/burgelab/miso/annotations/hg19_alt_events.zip
.. _3′READS annotations for (mm9) mouse genome (.zip): http://genes.mit.edu/burgelab/miso/annotations/tian_apa_events.zip
.. _Mus_musculus.NCBIM37.65.gff: http://genes.mit.edu/burgelab/miso/annotations/gene-models/Mus_musculus.NCBIM37.65.gff.zip
.. _Mus_musculus.NCBIM37.65.with_chr.gff: http://genes.mit.edu/burgelab/miso/annotations/gene-models/Mus_musculus.NCBIM37.65.with_chr.gff.zip
.. _Homo_sapiens.GRCh37.65.gff: http://genes.mit.edu/burgelab/miso/annotations/gene-models/Homo_sapiens.GRCh37.65.gff.zip
.. _Homo_sapiens.GRCh37.65.with_chr.gff: http://genes.mit.edu/burgelab/miso/annotations/gene-models/Homo_sapiens.GRCh37.65.with_chr.gff.zip
.. _Bowtie: http://bowtie-bio.sourceforge.net/
.. _Tophat: http://tophat.cbcb.umd.edu/
.. _IGV: http://www.broadinstitute.org/igv/
.. _Cufflinks: http://cufflinks.cbcb.umd.edu/
.. _PolyA DB: http://polya.umdnj.edu/polyadb/
.. _repository: https://github.com/yarden/MISO
.. _Perl script: http://seqanswers.com/forums/showthread.php?t=3201&highlight=GFF3
.. _miso-users: http://mailman.mit.edu/mailman/listinfo/miso-users
.. _White House adopts MISO: http://www.mediabistro.com/fishbowldc/white-house-soup-of-the-day-64_b53593#.Tp2c76k31tA.gmail
.. _bedtools: http://code.google.com/p/bedtools/
.. _sashimi_plot: sashimi.html

