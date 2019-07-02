API
===

This section documents the classes and functions of the libraries underlying ImpDAR. These really are the workhorses behind the executables that you would use for command-line processing. On the other hand, if you want to integrate the processing steps implemented by ImpDAR into another program, you will be interacting with these libraries.

The central component of ImpDAR processing is the :class:`~impdar.lib.RadarData.RadarData` class. Not only does this object store all the radar returns and auxiliary information, it also has a number of methods for processing.

Some processing steps may be implemented separately from the :class:`~impdar.lib.RadarData.RadarData` class. At present, just :func:`concatenation <impdar.lib.process.concat>`, is separate because it acts on multiple :class:`~impdar.lib.RadarData.RadarData` objects.

Contents:

.. toctree::
        :maxdepth: 2

        RadarData
        Plotting
        Picking
        load
        process
        ImpdarError
