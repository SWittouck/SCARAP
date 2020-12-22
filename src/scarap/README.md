# Architecture of scripts and functions

## Script architecture

Remark: all scripts import utils

* scarap.py:
    * commandline interface
    * imports taskwrappers
* taskwrappers.py:
    * wrappers for all tasks that check commandline arguments and dependencies
    * imports tasks_composite, tasks_core, checkers
* checkers.py:
    * functions that check commandline arguments and dependencies
* tasks_composite.py:
    * top-level code for composite tasks (tasks that depend on other tasks)
    * imports tasks_core
* tasks_core.py:
    * top-level code for core tasks (tasks that don't depend on other tasks)
    * imports readerswriters, computers, callers, pan
* pan.py
    * code for the various builtin pangenome inference methods
    * imports readerswriters, computers, callers
* computers.py:
    * functions that perform more complex computations
* readerswriters.py:
    * functions that perform simple reading/writing operations
* callers.py:
    * functions that call other software (e.g. hmmer tasks, mmseqs2 tasks, mafft)
* utils.py:
    * simple utility functions that are used by multiple scripts
    
## Architecture of pan.py

The script pan.py implements various pangenome inference strategies. Its functions are structured as follows:

* The top-level pangenome inference function is **infer_pangeome**. Its main arguments are a list of faa files and the pangenome inference strategy. It calls the function infer_superfamilies and then applies the function split_superfamily to all superfamilies in parallel. 
* The function **split_superfamily** initializes the necessary data structures depending on the requested strategy and calls a strategy-specific function with the name split_family_recursive_STR (where STR is the strategy name, e.g. split_family_recursive_LHT). 
* The **split_family_recursive_STR** function will split a superfamily in the set of final families, by splitting the family into two subfamilies recursively. Per recursion iteration, it goes through the following steps:
    1) it checks wether splitting if even an option (e.g. families with only one or two genes aren't splittable)
    2) it calls a strategy-specific splitting function to perform the actual split of the family into two subfamilies
    3) it checks if the split should actually happen based on the genome content of the subfamilies
    4) if the split should happen, it finishes the split and calls itself on the subfamilies to start the next iterations
* The functions that perform the actual split of a family into two subfamilies are called **split_family_STR**, where STR is again the name of the strategy. 

Some remarks about the implementation of the various strategies:

* The pangenome table ("pan" variable) is always passed to the recursive functions as a table that is indexed on the gene column. This makes it faster (I think) to subset it by genes, which is an operation that is frequently needed. 
* For consistency, all splitting functions (split_family_STR) return [pan1, pan2]. For some of them, it would be possible to just return [genes1, genes2], which would be simpeler and maybe even a bit faster (because the pan table technically only needs to be split if the splitting criterion tells us that the split will go through). However, I currently prefer to pay that small speed price to retain code consistency. 

## Family splitting strategies in pan.py

Four family splitting modules have been implemented in pan.py that can be combined to form a family splitting strategy:

* linclust (L): fixed-threshold clustering in linear time (mmseqs linclust)
* ficlin (F): clustering in a fixed number of clusters in linear time (own implementation using mmseqs align)
* hclust (H): multiple sequence alignment (MAFFT) followed by hierarchical clustering (scipy.cluster.hierarchy)
* tree (T): multiple sequence alignment (MAFFT) followed by phylogeny inference (IQTREE)

Three of the modules can be used to propose a binary split of a gene family: F, H and T; at least one of them needs to be part of the strategy. Three modules can be used to select representative sequences for the next module: L, F and H. Taking representatives can speed up the process and/or make it more scalable. An example of a strategy would be FH: selection of a fixed number of clusters by ficlin followed by hierarchical clustering of the representatives to determine the binary split of the family. 

By default, all strategies are lazy: they will attempt to re-use information from the parent, such as an hclust object or tree (however, for some strategies this is impossible). We indicate a non-lazy variant of a strategy with the suffix "-nl", e.g. H-nl (compute the full MSA and hierarchical clustering for each split). 
