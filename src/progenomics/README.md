# Architecture of scripts and functions

## Script architecture

Remark: all scripts import utils

* progenomics.py:
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
* readerswriters.py:
    * functions that perform simple reading/writing operations
* computers.py:
    * functions that perform more complex computations
* callers.py:
    * functions that call other software (e.g. hmmer tasks, mmseqs2 tasks, mafft)
* utils.py:
    * simple utility functions that are used by multiple scripts
    
## Architecture of pan.py

The script pan.py implements various pangenome inference strategies. Its functions are structured as follows:

* The top-level pangenome inference function is **infer_pangeome**. Its main arguments are a list of faa files and the pangenome inference strategy. It calls the function infer_superfamilies and then applies the function split_superfamily to all superfamilies in parallel. 
* The function **split_superfamily** initializes the necessary data structures depending on the requested strategy and calls a strategy-specific function with the name split_family_recursive_STR (where STR is the strategy name, e.g. split_family_recursive_LHT_F). 
* The **split_family_recursive_STR** function will split a superfamily in the set of final families, by splitting the family into two subfamilies recursively. Per recursion iteration, it goes through the following steps:
    1) it checks wether splitting if even an option (e.g. families with only one or two genes aren't splittable)
    2) it calls a strategy-specific splitting function to perform the actual split of the family into two subfamilies
    3) it checks if the split should actually happen based on the genome content of the subfamilies
    4) if the split should happen, it finishes the split and calls itself on the subfamilies to start the next iterations
* The functions that perform the actual split of a family into two subfamilies are called **split_family_STR**, where STR is again the name of the strategy. 

Some remarks about the implementation of the various strategies:

* The pangenome table ("pan" variable) is always passed to the recursive functions as a table that is indexed on the gene column. This makes it faster (I think) to subset it by genes, which is an operation that is frequently needed. Unfortunately, the purpose of this is a bit defeated in the linclust and hclust strategy components, because these frequently perform pd.merge operations, and those remove the index completely. Thus, the object is reindexed a number of times in each iteration. In an ideal world I would find a solution for this, but currently I don't think that the re-indexing is a speed bottleneck. 
* For consistency, all splitting functions (split_family_STR) return [pan1, pan2]. For some of them, it would be possible to just return [genes1, genes2], which would be simpeler and maybe even a bit faster (because the pan table technically only needs to be split if the splitting criterion tells us that the split will go through). However, I currently prefer to pay that small speed price to retain code consistency. 
