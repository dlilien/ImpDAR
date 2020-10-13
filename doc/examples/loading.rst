Loading Examples
================

There is not a lot documented here because loading supported files is extremely straightforward in ImpDAR. Loading (i.e. converting raw radar output into the ImpDAR/StoDeep matlab format) is accomplished in a single command with the :doc:`impdar load </../bin/impdar>` command.

The only real variation amongst filetypes is that you need to tell impdar what type of input file you are using. For example, for GSSI files, ``impdar load gssi fn [fn ...]`` will produce, for each input fn, a file with identical name but file ext '.mat'. Often, one might want to put all these outputs in a separate folder. The `-o folder_name` allows specification of an output folder. If you wanted to load PulseEkko data, it would be as simple as switching ``pe`` for ``gssi`` in the command above.
