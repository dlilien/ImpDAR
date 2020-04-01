.. ImpDAR documentation master file, created by
   sphinx-quickstart on Sun Jun  3 13:07:25 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ImpDAR's documentation!
==================================

ImpDAR is a flexible, open-source impulse radar processor that provides most of the benefits (and some additional features) compared to expensive commercial software. The starting point was the old St. Olaf deep radar matlab code. This code has a lot of history of contributors--I've tried to preserve acknowledgment of many of them in the file headers.

Support is gradually being added for a variety of file formats. Currently, GSSI, PulseEKKO, Radan, Blue Systems, DELORES, SEGY, gprMAX, seidart, Gecko, and legacy StoDeep files are supported. Available processing steps include various filtering operations, trivial modifications such as restacking, cropping, or reversing data, and a few different geolocation-related operations like interpolating to constant trace spacing. The primary interface is through the command line, which allows efficient processing of large volumes of data. An API, centered around the RadarData class, is also available to allow the user to use ImpDAR in other programs.

In addition to processing, ImpDAR can also be used for picking reflectors. Picking is generally an interactive process, and there is a light GUI for doing the picking. The GUI also provides support for basic processing operations, so you can see the effect of steps as you go along.

Requirements
------------

Python 2.7+ or 3.4+

`numpy <http://www.numpy.org>`_, 
`scipy <http://www.scipy.org>`_, 
`matplotlib <http://matplotlib.org>`_ 

To do anything involving geolocation, you will also need `GDAL <http://gdal.org>`_. The GUI, which is needed to be able to pick reflectors, requires `PyQt5 <https://pypi.org/project/PyQt5/>`_.
`SegYIO <https://github.com/equinor/segyio/>`_ is needed for SEGY support and for SeisUnix migration. `h5py <https://h5py.org>`_ is needed for some data formats.

.. include:: installation.rst

Contributing
------------
I would be thrilled to get pull requests for any additional functionality. In particular, it is difficult for me to add support for input formats for which I do not have example data--any development of readers for additional data types would be greatly appreciated.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
    
   installation.rst
   lib/index.rst
   bin/index.rst
   examples/index.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
