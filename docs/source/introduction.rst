Introduction
==============

Gimble is ...

Usage
***************

.. code-block:: bash

  usage: gimble <module> [<args>...] [-D -V -h]

  [Modules]
    preprocess            Preprocess input files
    setup                 Setup data store
    info                  Print information about DataStore
    blocks                Generate blocks from data in DataStore 
    windows               Generate windows from blocks in DataStore (requires blocks)
    query                 Query BED file of blocks (windows [TBI])
    model                 Build demographic model
    simulate              Simulate data [TBI] 
    makegrid              Make grid [TBI]
    gridsearch            Search grid [TBI]
    inference             Make inference [TBI] (requires blocks)
    
    partitioncds          Partition CDS sites in BED file by degeneracy in sample GTs 
    plotbed               Plot BED file [TBR]

.. literalinclude:: ../../cli/interface.py
   :start-after: """
   :end-before: """