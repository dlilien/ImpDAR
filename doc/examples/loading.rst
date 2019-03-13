Loading Examples
================

Loading (i.e. converting raw radar output into the ImpDAR/StoDeep matlab format) is accomplished in a single command with the :doc:`impdar load </../bin/impdar>`.

GSSI
----
``impdar load gssi fn [fn ...]`` will produce, for each input fn, a file with identical name but file ext '.mat'. Often, one might want to put all these outputs in a separate folder. The `-o folder_name` allows specification of an output folder.

PulseEkko
---------
``impdar load pe fn [fn ...]`` will produce, for each input fn, a file with identical name but file ext '.mat'. Often, one might want to put all these outputs in a separate folder. The `-o folder_name` allows specification of an output folder.
