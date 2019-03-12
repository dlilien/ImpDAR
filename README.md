# ImpDAR: an impulse radar processor

[![Build Status](https://travis-ci.org/dlilien/ImpDAR.svg?branch=master)](https://travis-ci.org/dlilien/ImpDAR) [![Build status](https://ci.appveyor.com/api/projects/status/uuef8aio2xbgiux8?svg=true)](https://ci.appveyor.com/project/dlilien/impdar)

This is a re-write of the St. Olaf Deep Radar processor in Python, adding some capability and pruning some dead limbs. This code has a lot of history of contributors--I've tried to preserve acknowledgment of many of them in the file headers. ImpDAR is intended to be more flexible than other available options. Support is gradually being added for a variety of file formats. Currently, GSSI and PulseEKKO files are supported. Available processing steps include various filtering operations, trivial modifications such as restacking, cropping, or reversing data, and a few different geolocation-related operations like interpolating to constant trace spacing.

In addition to processing, ImpDAR can also be used for picking reflectors. Picking is generally an interactive process, and there is something of a GUI for doing the picking. For the most up-to-date version of the picker, check out the "picker" branch, where most development of the picker is happening.

### Dependencies

#### Required
*Python 2 or 3* The package is tested on both [Windows](https://ci.appveyor.com/project/dlilien/impdar) and [Linux](https://travis-ci.org/dlilien/ImpDAR) on Python 2.7, 3.5, and 3.6. Other versions may work, but they are not tested.

[numpy](http://www.scipy.org)

[scipy](http://numpy.org)

[matplotlib](http://matplotlib.org)

#### Recommended
[GDAL](http://gdal.org) is needed to reproject out of WGS84, and thus for proper distance measurement.

[PyQt5](https://pypi.org/project/PyQt5/) is needed to run the GUI, which is needed for picking. You can do everything from the command line, and plot the results with matplotlib, without PyQt5. PyQt4 may work, but I haven't tested it.

## Documentation

Documentation of the various processing steps is [here](http://dlilien.github.io/ImpDAR). Examples are gradually being added to the documentation.

## Installation

Some explanation of other options is available in the main documentation, but I'll try to keep updating PyPi and so the easiest is `pip install impdar`.

## Contributing

I only have added support for different radar data types as needed--contributions for readers for other systems, whether commercial or custom, are always welcome.
