.. SCARAP documentation master file, created by
   sphinx-quickstart on Fri Oct 18 12:50:18 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

SCARAP: pangenome inference and comparative genomics of prokaryotes
====================================================================

.. image:: _static/img/scarap_logo.png
   :align: center
   :width: 200
   :alt: SCARAP logo

SCARAP is a toolkit with modules for various tasks related to comparative genomics of prokaryotes. SCARAP has been designed to be fast and scalable. Its main feature is pangenome inference, but it also has modules for direct core genome inference (without inferring the full pangenome), subsampling representatives from a (large) set of genomes and constructing a concatenated core gene alignment ("supermatrix") that can later be used for phylogeny inference. SCARAP has been designed for prokaryotes but should work for eukaryotic genomes as well. It can handle large genome datasets on a range of taxonomic levels; it has been tested on datasets with prokaryotic genomes from the species to the order level.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   getting-started
   scarap
   feedback
   citation
   license

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
