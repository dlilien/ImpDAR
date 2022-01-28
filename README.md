# ImpDAR: an impulse radar processor

[![DOI](https://zenodo.org/badge/134008583.svg)](https://zenodo.org/badge/latestdoi/134008583) [![Tests](https://github.com/dlilien/ImpDAR/actions/workflows/python-package.yml/badge.svg)](https://github.com/dlilien/ImpDAR/actions/workflows/python-package.yml) [![codecov](https://codecov.io/gh/dlilien/ImpDAR/branch/main/graph/badge.svg?token=R51QB61XNP)](https://codecov.io/gh/dlilien/ImpDAR)

ImpDAR is a suite of processing and interpretation tools for impulse radar (targeted for ice-penetrating radar applications but usable for ground-penetrating radar as well). The core processing steps and general terminology come from of the St. Olaf Deep Radar processor, but everything is re-written in python and significant improvements to speed, readability, documentation, and interface have been made across the board. However, this code has a lot of history of contributors--acknowledgment of many of them are preserved in the file headers.

ImpDAR is intended to be more flexible than other available options. Support is gradually being added for a variety of file formats. Currently, [GSSI](http://www.geophysical.com), [PulseEKKO](http://www.sensoft.ca), [Ramac](http://www.malagpr.com), [Blue Systems](http://www.bluesystem.ca/ice-penetrating-radar.html), DELORES, SEGY, [gprMAX](http://www.gprmax.com), Gecko, and legacy StoDeep files are supported. Available processing steps include various filtering operations, trivial modifications such as restacking, cropping, or reversing data, and a few different geolocation-related operations like interpolating to constant trace spacing. The integrated migration routines are in development but Stolt is working.

The primary interface to ImpDAR is through the command line, which allows efficient processing of large volumes of data. An API, centered around the RadarData class, is also available to allow the user to use ImpDAR in other programs.

In addition to processing, ImpDAR can also be used for interpreting the radargrams (i.e. picking layers). Picking is generally an interactive process, and there is a GUI for doing the picking; the GUI requires PyQt5, which may be annoying as a source build but is easy to install with Anaconda. The GUI also allows for basic processing steps, with the updates happening to the plot in real time. However, we have not added an 'undo' button to the GUI, so you are stuck with going back to your last saved version if you do not like the interactive results.

## Documentation

Documentation of the various processing steps is [here](https://impdar.readthedocs.io/en/latest/). There are examples of basic processing and plotting, and a longer example showing migration.

## Installation

Easiest is `pip install impdar`. Some explanation of other options is available in the main documentation, but PyPi will be updated with each release.

### Dependencies

#### Required
*Python 3* The package is tested on Python 3.7+. Older versions may work, but we have stopped testing on 2.7 since it has reached end of life. You can probably get 2.7 to work still, but no guarantees.

You also need:
[numpy](http://www.scipy.org),
[scipy](http://numpy.org),
[matplotlib](http://matplotlib.org),
[SegYIO](https://github.com/equinor/segyio/) for SEGY support and for SeisUnix migration, and
[h5py](https://h5py.org) is needed for some data formats.

I recommend just using Anaconda for your install, since it will also get you PyQt and therefore enable the GUI.

#### Recommended
[GDAL](http://gdal.org) is needed to reproject out of WGS84, and thus for proper distance measurement. Otherwise, distance calculations re going to moderately or severely incorrect.

[PyQt5](https://pypi.org/project/PyQt5/) is needed to run the GUI, which is needed for picking. You can do everything from the command line, and plot the results with matplotlib, without PyQt5.

Depending on whether you need migration routines, there may be some external dependencies. ImpDAR is designed to interface with [SeisUnix](http://https://github.com/JohnWStockwellJr/SeisUnix), which contains a number of powerful migration routines. You need to install SeisUnix yourself and get it on your path. If you are running windows, you need to figure out how to use Cygwin as well. However, the pure python migration routines in ImpDAR can work quite well, so don't let the difficulty of installing these compiled routines stop you from using those. ImpDAR searches for SeisUnix at the time of the call to the migration routine, so you can always add this later if you find that you need it.

## Contributing

Support for different radar data types has primarily been added as needed--contributions for data readers for other systems, whether commercial or custom, are always welcome.
