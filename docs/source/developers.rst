.. include:: <isogrk3.txt>

.. contents::

Developer guide to MISO 
=======================

This document is intended for people who want to modify, extend or alter the MISO source code. 
The focus will be on modifications to the core statistical inference engine.

Overview of MISO engine
=======================

Scoring functions. Single-end versus paired-end. Gibbs sampler, Metropolis Hastings, proposal
functions.


Coding guidelines
=================


Distribution guidelines
=======================

These guidelines are followed to simplify installation for users:

* **End-user installation must not depend on Cython.** The ``setup.py`` script is designed to compile the necessary
C/C++ source code (which could be custom code or Cython-generated), but it never runs Cython or deals with
``.pyx`` files. This is so that users don't have to install or use Cython in order to install MISO, and to prevent
situations where different Cython versions may produce different code. 
* **Cython-generated source files are not checked into repositories.*** This can create
situations where outdated files are distributed. The necessary C/C++ source code 
is compiled into the package using the ``MANIFEST.in`` file. This creates
source distributions that can be released that have the required C/C++ source files.
* **Numpy is not used from Cython source files.** All calls to ``numpy`` module are made on the 
pure Python-end, so that the Numpy C-API does not have to be used within Cython. Dependence on
the Numpy C-API can cause situations where Cython-generated code does not work with subsequent
Numpy releases, where the C-API might change. The Python Numpy API is less likely to change in an
incompatible way and is easier to modify if it does. Internally, the Cython source uses Cython or C arrays
rather than the Numpy ``np.ndarray`` array type.

The development testing/distribution pipeline is:

1. Create a source distribution by Cythonizing the necessary files. This step requires Cython, and is 
only intended for developers and not end-users: ::

   cd MISO/
   python setup.py sdist

This creates a source distribution that can be released or uploaded to Pypi. After this step, Cython is no
longer needed. 

2. Compile the source using ``python setup.py build_ext`` or using a package manager like ``pip`` ::

   pip install . 



