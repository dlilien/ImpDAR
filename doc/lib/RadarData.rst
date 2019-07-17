RadarData
=========
This page contains the documentation for the RadarData class, which is the basic object in ImpDAR.
If you are interacting with the API in a significant way, this is where you will find documentation from most of the things you care about, particularly how the data is stored and how to do basic processing steps on it.
All of the files to define the class are in impdar/lib/Radardata, with the basic initialization and class properties found in __init__.py and addional functionality spread across _RadarDataSaving, _RadarDataFiltering, and _RadarDataProcessing.

RadarData Base
--------------
.. autoclass:: impdar.lib.RadarData.RadarData
    :members: attrs_guaranteed, attrs_optional, chan, data, decday, dist, dt, lat, long, pressure, snum, tnum, trace_int, trace_num, travel_time, trig, trig_level, nmo_depth, elev, x_coord, y_coord, fn, check_attrs


Saving RadarData
----------------
These are all instance methods for saving information from a RadarData object.
They are defined in impdar/lib/RadarData/_RadarDataSaving.py.

.. automethod:: impdar.lib.RadarData.__init__.RadarData.save

.. automethod:: impdar.lib.RadarData.__init__.RadarData.save_as_segy

.. automethod:: impdar.lib.RadarData.__init__.RadarData.output_shp

.. automethod:: impdar.lib.RadarData.__init__.RadarData.output_csv


Processing RadarData
--------------------
These are all instance methods for processing data on a RadarData object.
They are defined in impdar/lib/RadarData/_RadarDataProcessing.py.

.. automethod:: impdar.lib.RadarData.__init__.RadarData.reverse

.. automethod:: impdar.lib.RadarData.__init__.RadarData.nmo

.. automethod:: impdar.lib.RadarData.__init__.RadarData.crop

.. automethod:: impdar.lib.RadarData.__init__.RadarData.restack

.. automethod:: impdar.lib.RadarData.__init__.RadarData.rangegain

.. automethod:: impdar.lib.RadarData.__init__.RadarData.agc

.. automethod:: impdar.lib.RadarData.__init__.RadarData.constant_space

.. automethod:: impdar.lib.RadarData.__init__.RadarData.elev_correct


Filtering Radar Data
--------------------
These are all instance methods for filtering data to remove noise.
They are defined in impdar/lib/RadarData/_RadarDataFiltering.py.

.. automethod:: impdar.lib.RadarData.__init__.RadarData.migrate

.. automethod:: impdar.lib.RadarData.__init__.RadarData.vertical_band_pass

.. automethod:: impdar.lib.RadarData.__init__.RadarData.adaptivehfilt

.. automethod:: impdar.lib.RadarData.__init__.RadarData.horizontalfilt

.. automethod:: impdar.lib.RadarData.__init__.RadarData.highpass

.. automethod:: impdar.lib.RadarData.__init__.RadarData.winavg_hfilt

.. automethod:: impdar.lib.RadarData.__init__.RadarData.hfilt

