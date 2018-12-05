Installation
============

 1. Install miniconda_.
 2. Install rdkit_ and dask_ using the following commands:
    
    - ``conda conda install -c rdkit rdkit``
    - ``conda install dask``
      
 3. Install the *insilico* library using the following command:
    
    - ``pip install git+https://github.com/SCM-NV/filterInsilico@master#egg=insilico-0.0.1``
  
Virtual Screening
=================
Virtual screening allows to filter a set of molecules based on their physico-chemical
properties from a potentially large number of candidates. 


.. _miniconda: http://conda.pydata.org/miniconda.html
.. _rdkit: http://www.rdkit.org
.. _dask: https://dask.org/
