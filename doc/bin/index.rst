Executables
===========

ImpDAR has four executables:

:doc:`impdar <impdar>` is a generic call that can process data, load data, or plot. Using this call, you can perform a number of processing steps in one go, saving time on loading and saving and saving disk space on not writing intermediate outputs.

:doc:`impproc <impproc>` is designed to give greater flexibility and cleaner syntax for processing. It only performs one processing step at a time, but will thus give you intermediate outputs, by default saved with names indicating the processing performed.

:doc:`impplot <impplot>` plots data, either as a radargram, as a line plot of power versus depth, or as the return power from a pick. It can either save the plot or bring it up for interactive panning and zooming.

:doc:`imppick <imppick>` calls up the interpretation GUI. Some processing can also be done through this GUI.

Contents:

.. toctree::
        :maxdepth: 2

        impdar
        impproc
        impplot
        imppick
