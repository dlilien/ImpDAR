ApresDiffData
=========
This page contains the documentation for the ApresDiffData class, which is the secondary object in ImpDAR's apres processor ``apdar``.
All of the files to define the class are in impdar/lib/ApresData, with the basic initialization and class properties found in __init__.py and addional functionality spread across _ApresDataSaving, _ApresDiffProcessing.

.. autoclass:: impdar.lib.ApresData.ApresDiffData
    :members: attrs_guaranteed, attrs_optional, data, data2, range, fn1, fn2, fn, unc1, unc2, ds, co, phi, w, w_err, w_0, eps_zz, bed, check_attrs

.. automethod:: impdar.lib.ApresData.ApresDiffData.save

.. automethod:: impdar.lib.ApresData.ApresDiffData.save

.. automethod:: impdar.lib.ApresData.ApresDiffData.phase_diff

.. automethod:: impdar.lib.ApresData.ApresDiffData.range_diff

.. automethod:: impdar.lib.ApresData.ApresDiffData.strain_rate

.. automethod:: impdar.lib.ApresData.ApresDiffData.bed_pick
