Installing SCARAP
=================

The easiest way to get started is to install SCARAP using conda.

Conda
-----

First, create and activate a new conda environment: ::

    conda create -n scarap python=3.11
    conda activate scarap

Then install from the recipe on bioconda: ::

    conda install bioconda::scarap    


Pip
---

First make sure that MAFFT and MMseqs2 are properly installed. Then install SCARAP with pip: ::

    pip install scarap


Manual install
--------------

You can also install SCARAP manually by cloning it and installing the following dependencies:

* `Python3 <https://www.python.org/>`_ version \>= 3.6.7 and < 3.13
    * `numpy <https://numpy.org/>`_ version \>= 1.16.5
    * `scipy <https://www.scipy.org/>`_ version \>= 1.4.1
    * `pandas version <https://pandas.pydata.org/>`_ \>= 1.5.3
    * `biopython <https://biopython.org/>`_ version \>= 1.67
    * `ete3 <http://etetoolkit.org/>`_ version \>= 3.1.1
* `MAFFT <https://mafft.cbrc.jp/alignment/software/>`_ version \>= 7.407
* `MMseqs2 <https://github.com/soedinglab/MMseqs2>`_ release 11, 12 or 13
