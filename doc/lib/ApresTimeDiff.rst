ApresTimeDiff
=========
This page contains the documentation for the Apres ApresTimeDiff class, which is an alternative object in ImpDAR's apres processor ``apdar``.
All of the files to define the class are in impdar/lib/ApresData, with the basic initialization and class properties found in __init__.py and addional functionality spread across _ApresDataSaving, _TimeDiffProcessing.

.. autoclass:: impdar.lib.ApresData.ApresTimeDiff
    :members: attrs_guaranteed, attrs_optional, data, data2, range, fn1, fn2, fn, unc1, unc2, ds, co, phi, w, w_err, w_0, eps_zz, bed, check_attrs

.. automethod:: impdar.lib.ApresData.ApresTimeDiff.save

.. automethod:: impdar.lib.ApresData.ApresTimeDiff.phase_diff

.. automethod:: impdar.lib.ApresData.ApresTimeDiff.range_diff

.. automethod:: impdar.lib.ApresData.ApresTimeDiff.strain_rate

.. automethod:: impdar.lib.ApresData.ApresTimeDiff.bed_pick
