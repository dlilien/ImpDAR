.. ImpDAR documentation master file, created by
   sphinx-quickstart on Sun Jun  3 13:07:25 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ImpDAR's documentation!
==================================

ImpDAR is a flexible, open-source impulse radar processor that provides most of the benefits (and some additional features) compared to expensive commercial software. The starting point was the old St. Olaf deep radar matlab code. This code has a lot of history of contributors--I've tried to preserve acknowledgment of many of them in the file headers.

Support is gradually being added for a variety of file formats. Currently, GSSI, PulseEKKO, Gecko, SEGY, and legacy StoDeep files are supported. Available processing steps include various filtering operations, trivial modifications such as restacking, cropping, or reversing data, and a few different geolocation-related operations like interpolating to constant trace spacing. The primary interface is through the command line, which allows efficient processing of large volumes of data. An API, centered around the RadarData class, is also available to allow the user to use ImpDAR in other programs.

In addition to processing, ImpDAR can also be used for picking reflectors. Picking is generally an interactive process, and there is something of a GUI for doing the picking. The GUI also provides support for basic processing operations, so you can see the effect of steps as you go along.

ImpDAR vs StODeep
-----------------

Compared to StoDeep, the jumping-off point for developing this software, ImpDAR provides all the same essential functionality but with a modern interface and version control. On the other hand, ImpDAR works on modern Python, requires no proprietary software, and has a user interface that should be easy to use for anybody who has familiarity with the command line and/or Python.
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

If you do not have a current (2.7 or 3+) python installation, you will need one to begin. I recommend `anaconda <https://anaconda.org/>`_. If you choose to use anaconda, you can install GDAL with `conda install -c conda-forge gdal`. GDAL is technically optional, but is needed for the majority of use cases.

If you want to be able to trace layers, you need PyQt5. PyQt5 ships with Anaconda, so with no further installation is needed if you installed python with Anaconda. Otherwise, you might need to build both Qt5 and PyQt5--you are on your own for those builds.

Now, you are ready to install impdar. You can get a version with

.. code-block:: bash

    pip install impdar
    
The downside here is that I need to update version tags in order to update on PyPi, and tagging releases is not a priority. To be sure that you have the newest version as a lot of development is happening, use

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

Contributing
------------
I would be thrilled to get pull requests for any additional functionality. In particular, it is difficult for me to add support for input formats for which I do not have example data--any development of readers for additional data types would be greatly appreciated.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
