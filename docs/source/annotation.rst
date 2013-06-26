.. include:: <isogrk3.txt>

.. contents::


Annotations for use with MISO and `sashimi_plot`_
=================================================

This page contains links to GFF annotations for use with MISO and `sashimi_plot`_. The GFF annotation format and how it is used by MISO is described in detail in the `MISO manual <index.html>`_.

Exon-centric annotations
------------------------

These annotations include GFF files (``.gff3`` extension) that can be used with MISO. 

Exon-centric annotations for human and mouse genomes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Version 1 of the human/mouse annotations (compiled 2008):
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* `Mouse genome (mm9) alternative events v1.0`_
* `Mouse genome (mm10) alternative events v1.0`_
* `Human genome (hg18) alternative events v1.0`_ 
* `Human genome (hg19) alternative events v1.0`_ 

These contain annotations of:

1. Skipped exons (SE)
2. Alternative 3'/5' splice sites (A3SS, A5SS)
3. Mutually exclusive exons (MXE)
4. Tandem 3' UTRs (TandemUTR)
5. Retained introns (RI)
6. Alternative first exons (AFE)
7. Alternative last exons (ALE)

Version 1 of the annotations for human and mouse genomes was derived from by Wang et. al. (2008) using ESTs and various annotation databases (like Ensembl, UCSC and AceView) to define alternative splicing events. Briefly, each splicing event was considered alternative if it was supported by several ESTs, and alternative tandem 3' UTRs (TandemUTR events) were derived from `PolyA DB`_.

Note that Version 1 of the annotations was originally made for mm9 and hg18, and the mm10 and hg19 annotation was made by coordinate mapping (using UCSC's ``liftOver`` utility) of mm9 to mm10, hg18 to hg19.

 .. warning::

    The lifted over Version 1 annotations of mm10/hg19 contain the ``ID`` entries in the GFF from mm9/hg18; however, the actual genomic coordinates, which are the only part read by MISO, have been lifted over to the more recent genomes. The ``ID`` value used in the GFF is arbitrary and is ignored by MISO; it is only used to encode the gene models hierarchy of genes, mRNAs and exons. Also note that lifting over is an imperfect process: not all events can always be fully lifted over.


**Mapping from alternative events to genes for Version 1 annotations**

Version 1 annotations from the links above contain a mapping from alternative events to genes, based on Ensembl annotation. These are tab-delimited files the first column (``event_id``) is the ``ID`` of the event from its GFF file and the second column (``gene_id``) corresponds to a comma-separated list of Ensembl identifiers for the gene(s) the event overlaps. If the event overlaps multiple genes (which could happen because multiple Ensembl identifiers are sometimes given to the same gene, or because the genes overlap and/or are contained within each other in the annotation), then multiple Ensembl identifiers will be listed. A mapping file is given for each event type (e.g. skipped exons, tandem 3' UTRs, etc.) Events that cannot be mapped to genes are recorded as ``NA``.


Version 2 (alpha release) of the human/mouse annotations (compiled June 2013):
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* `Mouse genome (mm9) alternative events v2.0`_
* `Mouse genome (mm10) alternative events v2.0`_
* `Human genome (hg18) alternative events v2.0`_
* `Human genome (hg19) alternative events v2.0`_

These contain annotations of:

1. Skipped exons (SE)
2. Alternative 3'/5' splice sites (A3SS, A5SS)
3. Mutually exclusive exons (MXE)
4. Retained introns (RI)

Version 2 of the annotations was derived by considering all transcripts annotated in Ensembl genes, knownGenes (UCSC) and RefSeq genes. The flanking exons to alternative exons were chosen using the "common shortest" rule, i.e. taking the shortest stretches of flanking that are most common among the annotated transcripts for the gene. The code used to generate these annotations is available as part of `rnaseqlib`_.

The annotations contain the following additional GFF attributes for each event's ``gene`` entry:

* ``ensg_id``: Ensembl ID for the gene the event falls within
* ``refseq_id``: RefSeq ID for the gene the event falls within
* ``gsymbol``: Gene symbol for the gene the event falls within

These annotations are still being tested. Comments on the annotation are welcomed. 


Exon-centric annotations for fly genome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* `Drosophila melanogaster alternative events (modENCODE)`_ (3.2M .zip)

These fly genome annotations were derived by the Graveley lab.


Isoform-centric annotations and reference gene models
-----------------------------------------------------

We provide GFF3 annotations based on UCSC Table Browser's version of Ensembl genes for the following genomes:

  * Mouse Ensembl genes from UCSC

    - Annotation for mm9: `mm9 ensGene GFF annotation`_

    - Annotation for mm10: `mm10 ensGene GFF annotation`_

  * Human Ensembl genes from UCSC

    - Annotation for hg18: `hg18 ensGene GFF annotation`_

    - Annotation for hg19: `hg19 ensGene GFF annotation`_

These can be used with MISO for isoform-centric quantitation, or with `sashimi_plot`_ to make plots of RNA-Seq data across gene models.

For convenience, we also provide GFF3 annotations of gene models from Ensembl (release 65), which were simply converted from Ensembl's GTF to GFF3 format and are otherwise identical to the Ensembl annotation.

  * Mouse Ensembl genes (16 M, .zip): `Mus_musculus.NCBIM37.65.gff`_
  * Human Ensembl genes (26 M, .zip): `Homo_sapiens.GRCh37.65.gff`_

Note that these annotations follow **Ensembl-style chromosome names** where as the UCSC-derived Ensembl annotations follow **UCSC-style chromosome names**.


Alternative 3' UTR annotations (hybrid)
---------------------------------------

In addition to exon-centric tandem 3' UTR annotations, alternative 3' UTR annotations for mouse (mm9) were made available by Wencheng Li and Bin Tian (these were derived by the 3' READS method: `Analysis of alternative cleavage and polyadenylation by 3′ region extraction and deep sequencing <http://www.nature.com/nmeth/journal/vaop/ncurrent/abs/nmeth.2288.html>`_). These contain two or more 3' UTR annotations per gene:

* `3′READS annotations for (mm9) mouse genome (.zip)`_


Updates
=======

**2013**:

* **Wed, Jun 26**: Released revised version of v1.0 annotations, where ALE event formatting error was fixed. Sanitized annotations to ensure start < end for hg19. In v2.0 annotations, an error in subset of RI event definitions was fixed.

.. _A primer on probabilistic inference: http://cocosci.berkeley.edu/tom/papers/tutorial.pdf
.. _rnaseqlib: http://rnaseqlib.readthedocs.org
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
.. _Mouse genome (mm9) alternative events v1.0: http://genes.mit.edu/burgelab/miso/annotations/miso_annotations_mm9_v1.zip
.. _Mouse genome (mm10) alternative events v1.0: http://genes.mit.edu/burgelab/miso/annotations/miso_annotations_mm10_v1.zip
.. _Human genome (hg18) alternative events v1.0: http://genes.mit.edu/burgelab/miso/annotations/miso_annotations_hg18_v1.zip
.. _Human genome (hg19) alternative events v1.0: http://genes.mit.edu/burgelab/miso/annotations/miso_annotations_hg19_v1.zip
.. _Mouse genome (mm9) alternative events v2.0: http://genes.mit.edu/burgelab/miso/annotations/ver2/miso_annotations_mm9_v2.zip
.. _Mouse genome (mm10) alternative events v2.0: http://genes.mit.edu/burgelab/miso/annotations/ver2/miso_annotations_mm10_v2.zip
.. _Human genome (hg18) alternative events v2.0: http://genes.mit.edu/burgelab/miso/annotations/ver2/miso_annotations_hg18_v2.zip
.. _Human genome (hg19) alternative events v2.0: http://genes.mit.edu/burgelab/miso/annotations/ver2/miso_annotations_hg19_v2.zip
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



