ApresData
=========
This page contains the documentation for the ApresData class, which is the basic object in ImpDAR's apres processor ``apdar``.
All of the files to define the class are in impdar/lib/ApresData, with the basic initialization and class properties found in __init__.py and addional functionality spread across _ApresDataSaving, _ApresDataProcessing.

ApresData Base
--------------
.. autoclass:: impdar.lib.ApresData.ApresData
    :members: attrs_guaranteed, attrs_optional, data, decday, dt, snum, cnum, bnum, chirp_num, chirp_att, chirp_time, travel_time, frequencies, lat, long, x_coord, y_coord, elev temperature1, temperature2, battery_voltage, Rcoarse, uncertainty, fn, check_attrs


Saving ApresData
----------------
These are all instance methods for saving information from a ApresData object.
They are defined in impdar/lib/ApresData/_ApresDataSaving.py.

.. automethod:: impdar.lib.ApresData.__init__.ApresData.save

Processing ApresData
--------------------
These are all instance methods for processing data on a ApresData object.
They are defined in impdar/lib/ApresData/_ApresDataProcessing.py.

.. automethod:: impdar.lib.ApresData.__init__.ApresData.apres_range

.. automethod:: impdar.lib.ApresData.__init__.ApresData.phase_uncertainty

.. automethod:: impdar.lib.ApresData.__init__.ApresData.stacking
