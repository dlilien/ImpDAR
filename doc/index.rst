.. ImpDAR documentation master file, created by
   sphinx-quickstart on Sun Jun  3 13:07:25 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ImpDAR's documentation!
==================================

ImpDAR is a radar processor based on the St Olaf Deep radar processor. It allows for filtering, plotting, distance- and depth-correcting on a variety of file formats.  It does not yet have the full functionality of the original, but the most common steps are incorporated.

ImpDAR vs StODeep
-----------------

Currently, ImpDAR incorporates all major elements of StoDeep, but none of StoInterpret. On the other hand, ImpDAR works on modern Python, requires no proprietary software, and has a user interface that should be easy to use for anybody who has familiarity with the command line and/or Python. I would be thrilled to get pull requests for remaining missing functionality (particularly additional input formats).

Requirements
------------

Python 2 or 3, 
`numpy <http://www.numpy.org>`_, 
`scipy <http://www.scipy.org>`_, 
`matplotlib <http://matplotlib.org>`_ 

To do anything involving geolocation, you will also need `GDAL <http://gdal.org>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   lib/index.rst
   bin/index.rst

Installation
------------

If you do not have a current (2.7 or 3+) python installation, you will need one to begin. I recommend `anaconda <https://anaconda.org/>`_, though for users of homebrew a current homebrew python may be more convenient. Anaconda comes with the scientific python stack, but with homebrew you will need to `brew install numpy scipy; pip install matplotlib`.

If you choose to use anaconda, you can install GDAL with `conda install -c conda-forge gdal`. If you are using homewbrew, `brew install gdal` should do the trick. GDAL is technically optional, but is needed for the majority of use cases.

Finally, you are ready to install impdar. You can get a version with

.. code-block:: bash

    pip install impdar
    
The only downside here is that I need to remember to update version tags in order to update on PyPi, and I may not remember that often. To be sure that you have the newest version, use

.. code-block:: bash

    git clone https://github.com/dlilien/ImpDAR.git
    cd impdar
    python setup.py install

This requires `git <https://git-scm.com/downloads>`_. 

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
