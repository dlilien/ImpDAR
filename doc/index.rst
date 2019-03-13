.. ImpDAR documentation master file, created by
   sphinx-quickstart on Sun Jun  3 13:07:25 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ImpDAR's documentation!
==================================

This is a re-write of the St. Olaf Deep Radar processor in Python, adding some capability and pruning some dead limbs. This code has a lot of history of contributors--I've tried to preserve acknowledgment of many of them in the file headers. ImpDAR is intended to be more flexible than other available options. Support is gradually being added for a variety of file formats. Currently, GSSI and PulseEKKO files are supported. Available processing steps include various filtering operations, trivial modifications such as restacking, cropping, or reversing data, and a few different geolocation-related operations like interpolating to constant trace spacing. The primary interface is through the command line, which allows efficient processing of large volumes of data. An API, centered around the RadarData class, is also available to allow the user to use ImpDAR in other programs.

In addition to processing, ImpDAR can also be used for picking reflectors. Picking is generally an interactive process, and there is something of a GUI for doing the picking. For the most up-to-date version of the picker, check out the "picker" branch, where most development of the picker is happening.

ImpDAR vs StODeep
-----------------

Currently, ImpDAR incorporates all major elements of StoDeep except migration, and the essential pieces of StoInterpret. On the other hand, ImpDAR works on modern Python, requires no proprietary software, and has a user interface that should be easy to use for anybody who has familiarity with the command line and/or Python. I would be thrilled to get pull requests for remaining missing functionality (particularly additional input formats).

Requirements
------------

Python 2.7+ or 3.4+, 
`numpy <http://www.numpy.org>`_, 
`scipy <http://www.scipy.org>`_, 
`matplotlib <http://matplotlib.org>`_ 

To do anything involving geolocation, you will also need `GDAL <http://gdal.org>`_. The GUI, which is needed to be able to pick reflectors, requires `PyQt5 <https://pypi.org/project/PyQt5/>`_.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   lib/index.rst
   bin/index.rst
   examples/index.rst

Installation
------------

If you do not have a current (2.7 or 3+) python installation, you will need one to begin. I recommend `anaconda <https://anaconda.org/>`_, though for users of homebrew a current homebrew python may be more convenient. Anaconda comes with the scientific python stack, but with homebrew you will need to `brew install numpy scipy; pip install matplotlib`.

If you choose to use anaconda, you can install GDAL with `conda install -c conda-forge gdal`. If you are using homewbrew, `brew install gdal` should do the trick. GDAL is technically optional, but is needed for the majority of use cases.

If you want to be able to trace layers, you need PyQt5. PyQt5 ships with Anaconda, so with no further installation is needed if you installed python with Anaconda. Otherwise, you might need to build both Qt5 and PyQt5.

Finally, you are ready to install impdar. You can get a version with

.. code-block:: bash

    pip install impdar
    
The downside here is that I need to update version tags in order to update on PyPi, and that may not happen often. To be sure that you have the newest version as a lot of development is happening, use

.. code-block:: bash

    git clone https://github.com/dlilien/ImpDAR.git
    cd impdar
    python setup.py install

This requires `git <https://git-scm.com/downloads>`_. 

If you want to run the bleeding-edge version with the picker, use 
.. code-block:: bash

    git clone https://github.com/dlilien/ImpDAR.git
    git checkout picker
    cd impdar
    python setup.py install

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
