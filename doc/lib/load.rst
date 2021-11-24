Loading data
============

These are functions for loading loading radar data, generally from raw formats, to be used in a program or saved in ImpDAR's .mat format and used later.

For every filetype that ImpDAR can handle (e.g. GSSI .DZT files, gprMax .h5 files), there is a dedicated file for loading that filetype in `impdar/lib/load`. These files generally define a single method, which returns an `impdar.lib.RadarData.RadarData` object, with information specific to the filetype loaded in. The user does not need to interact with these files (unless they need to add functionality). For some of the systems, documentation is sparse and this is a challenge, while for others documentation is readily available (e.g. `Blue Systems <https://iprdoc.readthedocs.io/en/latest/>`).

Instead, to load data for interactive use, a generic `load` command, which takes a filetype as an argument, is defined in `impdar.lib.load.__init__`. This wrapper provides some conveniences for handling multiple files as well. There is also a `load_and_exit` command in that file, which can be used if the user does not want to interact with the data at load time, but wants the filetype converted to ImpDAR's .mat for convenience.

.. automethod:: impdar.lib.load.load

.. automethod:: impdar.lib.load.load_and_exit
