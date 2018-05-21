# impdar: an impulse radar processor

This is a re-write of the St. Olaf Deep Radar processor in Python, adding some capability and pruning some dead limbs. This code has a lot of history of contributors: I've tried to preserve acknowledgment of many of them in the file headers.


## Installation

### Dependencies

#### Required
**Python 2 or 3** Sticking with 3.6 or greater is prudent
**SCIPY, NUMPY, and MATPLOTLIB**

#### Recommended
**GDAL** Needed to reproject out of WGS84, and thus for proper distance measurement
