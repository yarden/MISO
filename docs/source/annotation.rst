.. include:: <isogrk3.txt>

.. contents::


Annotations for use with MISO and `sashimi_plot`_
=================================================

This page contains links to GFF annotations for use with MISO and `sashimi_plot`. The annotations
are described in detail in the `MISO manual <index.html>`_.

Exon-centric
------------

Available for human/mouse/fly genomes:

* `Drosophila melanogaster alternative events (modENCODE)`_ (3.2M .zip)
* `Mouse genome (mm9) alternative events`_ (7.9M, .zip)
* `Human genome (hg18) alternative events`_ (45M, .zip)
* `Human genome (hg19) alternative events`_ (8.6M, .zip)


Gene and isoform-centric annotations
------------------------------------

We provide GFF3 annotations based on UCSC Table Browser's version of Ensembl genes for the following genomes:

  * Mouse Ensembl genes from UCSC

    - Annotation for mm9: `mm9 ensGene GFF annotation`_

    - Annotation for mm10: `mm10 ensGene GFF annotation`_

  * Human Ensembl genes from UCSC

    - Annotation for hg18: `hg18 ensGene GFF annotation`_

    - Annotation for hg19: `hg19 ensGene GFF annotation`_

These can be used with MISO for isoform-centric quantitation, or with `sashimi_plot` to make plots of RNA-Seq data across gene models.

For convenience, we also provide GFF3 annotations of gene models from Ensembl (release 65), which were simply converted from Ensembl's GTF to GFF3 format and are otherwise identical to the Ensembl annotation.

  * Mouse Ensembl genes (16 M, .zip): `Mus_musculus.NCBIM37.65.gff`_
  * Human Ensembl genes (26 M, .zip): `Homo_sapiens.GRCh37.65.gff`_

Note that these annotations follow **Ensembl-style chromosome names** where as the UCSC-derived Ensembl annotations follow **UCSC-style chromosome names**.


Alternative 3' UTR annotations (hybrid)
---------------------------------------

In addition to exon-centric tandem 3' UTR annotations, alternative 3' UTR annotations for mouse (mm9) were made available by Wencheng Li and Bin Tian (these were derived by the 3' READS method). These contain two or more 3' UTR annotations per gene:

* `3′READS annotations for (mm9) mouse genome (.zip)`_

Reference gene annotations
--------------------------

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



