impproc
=======

An executable to perform single processing steps.

This has a lot of convenience in terms of the call since you get more help with commands, more control of arguments, control over the order in which things are done, etc, but has the disadvantage of requiring a call/load/write for every step.

You can get a list of commands with ``impproc -h``

For any individual command, you can get more help by running ``impproc [command] -h``.

Examples
--------

A sample workflow might be something like


.. code-block:: bash

    # make directories for the output
    mkdir bandpass hfilt nmo

    # Vertical bandpass from 150-450MHz (loading in the raw data with the -gssi flag)
    impproc vbp 150 450 -gssi *.DZT -o bandpass/

    # do some horizontal filtering on that output
    impproc hfilt 1000 2000 bandpass/*.mat -o hfilt

    # finally do a conversion to the time domain
    impproc nmo 10 hfilt/*.mat -o nmo
    

The same processing steps can be done without separating the output into different folders. At risk of file clutter, the workflow could be

.. code-block:: bash

    # Vertical bandpass from 150-450MHz (loading in the raw data with the -gssi flag)
    impproc vbp 150 450 -gssi *.DZT

    # do some horizontal filtering on that output
    impproc hfilt 1000 2000 *_vbp.mat

    # finally do a conversion to the time domain
    impproc nmo 10 *_hfilt.mat

    # Outputs are now sitting around with _vbp_hfilt_nmo before the extension

A similar example, with visualization of the outputs, is :doc:`here </../examples/processing>`.

Usage
-----

.. argparse::
    :module: impdar.bin.impproc
    :func: _get_args
    :prog: impproc
