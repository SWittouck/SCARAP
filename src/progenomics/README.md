# Scripts

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
    * imports readerswriters, computers, callers
* readerswriters.py:
    * functions that perform simple reading/writing operations
* computers.py:
    * functions that perform more complex computations
* callers.py:
    * functions that call other software (e.g. hmmer tasks, mmseqs2 tasks, mafft)
* utils.py:
    * simple utility functions that are used by multiple scripts