# impdar: an impulse radar processor

[![Build Status](https://travis-ci.org/dlilien/impdar.svg?branch=master)](https://travis-ci.org/dlilien/impdar)

This is a re-write of the St. Olaf Deep Radar processor in Python, adding some capability and pruning some dead limbs. This code has a lot of history of contributors: I've tried to preserve acknowledgment of many of them in the file headers.

### Dependencies

#### Required
*Python 2 or 3* Sticking with 3.6 or greater is prudent

[numpy](http://www.scipy.org),
[scipy](http://numpy.org),
[matplotlib](http://matplotlib.org)

#### Recommended
*GDAL* Needed to reproject out of WGS84, and thus for proper distance measurement

## Documentation

Reasonably complete documentation is [here](http://dlilien.github.io/impdar)
