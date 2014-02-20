.. include:: <isogrk3.txt>

Updates
=======

**2014**

* Renamed scripts. The core MISO scripts are now named ``miso``, ``summarize_miso`` and ``compare_samples``. Python extensions (``.py``) have been dropped. Entry points are used in Python installations to resolve issues with scripts not being registered as executables after installation.

* Added ``--no-wait`` option to ``miso`` to prevent main process from waiting on cluster jobs

* Fixed bug that prevented ``--apply-both`` from working in ``filter_events``

**2013**

* **Jun 26**: Fixed formatting error for version 1 ALE events. Reorganized documentation.

* **Apr 26**: Released version ``0.4.9``. Contains Sashimi plot-related bug fixes.

* **Apr 14**: Released version ``0.4.8``. Includes several important bug fixes and a new much faster version of the inference algorithm for single-end reads:

  * New, much faster inference scheme for single-end reads (implemented by Gabor Csardi) 

  * Fixed bug in parsing of Ensembl style chromosomes

  * Fixed bug in ``--prefilter`` option. Added ``-p`` option to ``run_events_analysis.py`` for determining the number of processors to be used.

  * Fixed potential deadlock situation when running locally on machines with multiple cores.

  * Added ``miso_zip`` utility for compressing/uncompressing MISO output.

  * Better error checking for incompatibilities between BAM and GFF files, e.g. cases where the chromosome naming conventions between the annotation and the BAM differ.

  * Fixed entries in hg19 alternative event annotations where start > end. 

  * Changed default MISO settings file to omit queue names from cluster submissions.

* **Tue, Jan 1**: Released version ``0.4.7``. Includes several new features, including:

  * Support for strand-specific reads
  
  * Support for running MISO locally using multi-cores

  * Prefilter feature: prefilter events with low coverage on long runs (requires Bedtools)


**2012**

* **Thu, Dec 20**: Posted new mouse genome annotation for alternative cleavage and polyadenylation events (TandemUTR) and alternative last exon events (ALE) from `Analysis of alternative cleavage and polyadenylation by 3′ region extraction and deep sequencing <http://www.nature.com/nmeth/journal/vaop/ncurrent/abs/nmeth.2288.html>`_ (Nature Methods, Dec. 2012), see :ref:`tian-apa`. Thanks to Wencheng Li and Bin Tian for creating this annotation!

* **Tue, Dec 4**: Posted mappings from alternative events to genes (see :ref:`events-to-genes`) and GFF annotations for isoform-centric inference (see :ref:`iso-centric`).

* **Thu, Sept 27**: Released ``misopy-0.4.6``. This release fixes a packaging error with ``sashimi_plot`` (test case files were omitted.) No changes were made to MISO. Thanks to Schragi Schwartz and Rahul Satija.

* **Tue, Sept 4**: Released ``misopy-0.4.5``. Fixed error in parsing of settings file.

* **Thur, Jul 26**: Released ``misopy-0.4.4``. Turned off autoconvergence features of sampler that caused delayed run times on certain events. (Also introduced a parameter ``num_chains`` in MISO settings file to control the number of independent MCMC chains used by the sampler.)

* **Fri, Jun 1**: Released ``misopy-0.4.3``. Included one main bug fix to an instability issue related to Bayes factors which affected low coverage events and better error handling for ``sashimi_plot``. The Bayes factor bug resulted in some very low coverage events having highly variable Bayes factor values across MISO runs. This issue does not affect the |Psi| values or their confidence intervals, only the Bayes factor values. Events that meet most standard coverage filters applied to read counts will not be affected. Thanks to Ersen Kavak, Essi Laajala and Schragi Schwartz. 

* **Fri, May 4**: Release of ``misopy-0.4.2``. Included some bug fixes, notably:

  - Showing complete output for events 
  - Fixed bug that prevented paired-end end inserts where the read mates overlap from being aligned 
  - Proper parsing of event names with multiple periods in their GFF ``ID`` field
  - Made SGE cluster submission scripts automatically executable 

  Thanks to Rob Bradley, Eric Suh, Brad Friedman and others for their comments.

* **Wed, Feb 1**: Release of ``misopy-0.4.1``. This is mainly a feature release. Updates include gene information being outputted in summary files and sample comparison files, dynamic scaling for ``sashimi_plot`` (described in its documentation page), options for labeling summary/comparison files, and other minor features. Bugfixes include: correct inclusion of test files with ``pypi`` releases, correction of edge cases that caused ``compare-samples`` to fail in case a multi-isoform event was improperly sampled.

* **Fri, Jan 13**: Release ``misopy-0.4``. Includes mainly fix to serious bug that caused MISO to drop junction reads for events containing alignments with insertions (e.g. CIGAR strings containing ``I``). Users of Tophat and other read mappers that allow insertions in alignments please take notice. 

* **Sun, Jan 8**: Release ``misopy-0.3``. Includes mainly ``sashimi_plot`` features/bugfixes.

* **Tue, Jan 3**: Released ``misopy-0.2``. Includes bugfixes to test cases and new ``sashimi_plot`` features.

**2011**

* **Mon, Dec 26**: We've done some major reorganization of MISO and ``sashimi_plot`` to make them official Python packages. MISO can now be installed like as a standard Python package from the Python Package Index (see our `pypi`_ page.) The manual has been changed to reflect the new and improved installation procedure, which makes it easier to invoke MISO from the command line, requiring no modifications to your environment path variables.  

    .. note::
      The new reorganized version of MISO cannot be used with older indexed/pickled files, and requires GFF files to be **re-indexed** by calling ``index_gff.py`` from the new package. Indexed/pickled files created with older versions of MISO will not work with the new version. We apologize for the inconvenience. Future updates to MISO will not need this re-indexing step, but our refactoring of the code results in this being required for this upgrade.


* **Fri, Dec 16**: Added support for cluster submission using SGE -- thanks to Michael Lovci for implementing this feature. Fixed edge case that caused ``filter_events`` to fail, thanks to Kuan-Ting Lin for pointing this out.

* **Wed, Dec 14**: Several new updates, including:

  - Utilities for computing and plotting insert-length distributions. This makes use of the awesome `bedtools`_ software written by Aaron Quinlan. Special thanks to Aaron for implementing our feature request, now available as the ``-intervals`` feature in ``tagBam``, which we make heavy use of. Users of `bedtools`_ should check it out!

  - Relaxed restrictions on GFF parser: can now accept "transcript" in addition to "mRNA" entries, and gracefully skip over genes with no mRNAs/transcripts in the file.

  - Posted GFF3 versions of Ensembl gene models for human/mouse to be used with isoform-centric analyses.

* **Sun, Dec 11**: Removed platform-specific pickled files from annotation zip files.

* **Sat, Dec 3**: `sashimi_plot <sashimi.html>`_, a tool for visualizing RNA-Seq reads along exons/junctions as well as MISO output is released! 

* **Fri, Dec 2**: Underscores (``_``) are now allowed in GFF annotations. 

* **Fri, Nov 11 (11/11/11)**: We posted a link to MISO annotations of alternative events for the Drosophila melanogaster genome, courtesy of Brent Graveley and the modENCODE project. Thanks to Brent for compiling and sharing these helpful annotations! We encourage those who work on Drosophila transcriptome datasets to try these and contribute more annotations of novel alternative RNA processing events.

* **Thu, Oct 27**: Released ``fastmiso``, a version of MISO written in the C programming language which is 60-100x faster than the Python only version.  The MISO version contains a Python interface that is identical to the original Python-onlyMISO version, and thus will require no modification on the part of users aside from compiling the C code. An interface to MISO in the R language is also available (currently undocumented.) 

* **Tue, Oct 18**: The `White House adopts MISO`_. From the article: *"Miso would be so upset if I didn’t tell you what it was," said Host Chuck Todd employing one of his many puns.*

* **Fri, Sept 9**: Human alternative event annotations (in GFF format) for hg19 are posted. Thanks to Brent Graveley for generating these.

* **Sun, Sept 4**: Fixed bug in GFF parser (``GFF.py``) that incorrectly skipped gene IDs containing DOS /  Windows newline character (``\r``).  Thanks to Liu Yujing for pointing this out.

* **Mon, Aug 29**: On some filesystems, the ``index_gff.py`` script outputted files with names too long for the operating system when indexing the ``ALE`` gff events (alternative last exons.) A new set of ALE ``.gff3`` files and their indexed versions was uploaded that avoids this problem. This should only affect ALE events. Thanks to Chia-Ho Lin for pointing this out.

* **Mon, Aug 15**: Formatting error in "ID" field of AFE/ALE events was corrected. New ``.gff3`` files and their pickled versions were posted. Thanks to Jessica Hurt.

* **Tues, Jul 12**: Bug in paired-end computation has been corrected. Thanks to Sol Katzman for documenting the bug and coming up with a test case. Also, mailing list for MISO users is now available: `miso-users`_ (http://mailman.mit.edu/mailman/listinfo/miso-users)

* **Mon, Jul 11**: Bug in alignment of BAM reads to certain gene models has been fixed (thanks to Gábor Csárdi for patch.) This bug might have preferentially affected "retained intron" alternative splicing events. Paired-end reads now support overhang constraint (which will be applied to both read mates.)

* **Sat, Mar 26**: Fixed bug in paired-end assignment of reads. Changes from ``dev`` branch are now merged into the ``master`` branch of the GitHub repository. Thanks again to Marvin Jens and Eric T. Wang.

* **Tue, Mar 8**: Updated human/mouse splicing event annotation zip files to contain both raw gff3 files as well as their indexed (pickled) versions. Thanks to Marvin Jens.

* **Wed, Feb 23**: The MISO GitHub repository has been relocated to this `repository`_ (``git@github.com:yarden/MISO.git``).

* **Tue, Feb 22**: Incorporated Rory Kirchner's changes to ``filter_events.py`` (from his fork of the MISO repository.) Supports filtering on biological replicates with the ``--votes`` option. Thanks to Rory for his code!

* **Fri, Feb 18**: Added warning to detect cases where sampling might not have been successful or might not have mixed.

* **Wed, Feb 9**: Added ``filter_events.py``, which allows easy filtering
  of events based on coverage and magnitude of change from a MISO Bayes
  factor comparison file.

**2010**

* **Thurs, Dec 30:** Implemented faster indexing scheme of GFF genes,
  which should significantly reduce MISO run times when using GFF
  annotations. To use the new indexing scheme, the **index must be
  rebuilt** using ``index_gff.py``. Also added overhang support for (single-end) BAM
  files. These additions are available in the "dev" branch of the GitHub repository.

* **Wed, Dec 22:** Fixed chromosome naming bug in SAM parsing.

* **Tues, Dec 21:** Critical filename parsing bug in ``--compare-samples``
  option fixed. Please rerun ``--compare-samples`` on MISO results to
  obtain corrected output. Much thanks to Sol Katzman (UCSC) for catching the bug and
  offering a solution.

* **Mon, Dec 20:** MISO now also accepts GFF files that do not set the
  "ID" field of exon entries. These entries get assigned a default
  ID. Edge case related to short isoforms also fixed. Completed
  support for Bayes factors and delta posterior densities for genes
  with arbitrarily many isoforms. The Github `repository`_ for MISO is
  now public

* **Thurs, Dec 2:** Fixed isoform length edge case. Removed unnecessary
  plotting-related code.

* **Mon, Nov 15:** Fixed bug that caused crash for GFF annotations that
  give genes with only a single isoform (only one "mRNA" entry). These
  events/genes are now skipped. Added explicit error logging. 

* **Sun, Nov 14:** Updated ``index_gff`` script so that GFF
  annotations of mRNAs are indexed in the order in which they appear
  in the file. This means that if the first mRNA in the file is A, second is B,
  and third is C, then |Psi| will be a 3-element vector with the
  probabilities of A, B and C respectively. New versions of pickled
  GFF event files were uploaded (for mm9 and hg18) to reflect this.
