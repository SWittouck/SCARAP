Installing SCARAP
=================

You can install SCARAP through `conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`_ ::
    
    git clone https://github.com/swittouck/scarap.git
    cd scarap
    conda env create -f environment.yml 

You can then run SCARAP as follows ::
    
    conda activate scarap
    scarap -h
    conda deactivate

You can also install SCARAP manually by cloning it and installing the following dependencies:

* `Python3 <https://www.python.org/>`_ version \>= 3.6.7
    * `numpy <https://numpy.org/>`_ version \>= 1.16.5
    * `scipy <https://www.scipy.org/>`_ version \>= 1.4.1
    * `pandas version <https://pandas.pydata.org/>`_ \>= 1.5.3
    * `biopython <https://biopython.org/>`_ version \>= 1.67
    * `ete3 <http://etetoolkit.org/>`_ version \>= 3.1.1
* `MAFFT <https://mafft.cbrc.jp/alignment/software/>`_ version \>= 7.407
* `MMseqs2 <https://github.com/soedinglab/MMseqs2>`_ release 11, 12 or 13