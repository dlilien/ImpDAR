ApresDiff
=========
This page contains the documentation for the ApresData class, which is the basic object in ImpDAR's apres processor ``apdar``.
All of the files to define the class are in impdar/lib/ApresData, with the basic initialization and class properties found in __init__.py and addional functionality spread across _ApresDataSaving, _ApresDataProcessing.

ApresDiff Base
--------------
.. autoclass:: impdar.lib.ApresData.ApresDiff
    :members: attrs_guaranteed, attrs_optional, data, data2, Rcoarse, fn1, fn2, fn, unc1, unc2, ds, co, phi, w, w_err, w_0, eps_zz, bed, check_attrs


Saving ApresDiff
----------------
These are all instance methods for saving information from a ApresDiff object.
They are defined in impdar/lib/ApresData/_ApresDataSaving.py.

.. automethod:: impdar.lib.ApresData._ApresDataDifferencing.ApresDiff.save

Processing ApresDiff
--------------------
These are all instance methods for processing data on a ApresDiff object.
They are defined in impdar/lib/ApresData/_ApresDataDifferencing.py.

.. automethod:: impdar.lib.ApresData._ApresDataDifferencing.ApresDiff.phase_diff

.. automethod:: impdar.lib.ApresData._ApresDataDifferencing.ApresDiff.range_diff

.. automethod:: impdar.lib.ApresData._ApresDataDifferencing.ApresDiff.strain_rate

.. automethod:: impdar.lib.ApresData._ApresDataDifferencing.ApresDiff.bed_pick
