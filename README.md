# impdar: an impulse radar processor

[![Build Status](https://travis-ci.org/dlilien/ImpDAR.svg?branch=master)](https://travis-ci.org/dlilien/ImpDAR) [![Build status](https://ci.appveyor.com/api/projects/status/uuef8aio2xbgiux8?svg=true)](https://ci.appveyor.com/project/dlilien/impdar)

This is a re-write of the St. Olaf Deep Radar processor in Python, adding some capability and pruning some dead limbs. This code has a lot of history of contributors: I've tried to preserve acknowledgment of many of them in the file headers.

### Dependencies

#### Required
*Python 2 or 3* The package is tested on both [Windows](https://ci.appveyor.com/project/dlilien/impdar) and [Linux](https://travis-ci.org/dlilien/ImpDAR) on Python 2.7, 3.5, and 3.6. Other versions may work, but they are not tested.

[numpy](http://www.scipy.org),
[scipy](http://numpy.org),
[matplotlib](http://matplotlib.org)

#### Recommended
*GDAL* Needed to reproject out of WGS84, and thus for proper distance measurement.

## Documentation

Reasonably complete documentation is [here](http://dlilien.github.io/ImpDAR)
