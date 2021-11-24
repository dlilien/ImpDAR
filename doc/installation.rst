Installation
------------

Beginner
________

If you do not have a current (3.6+) python installation, you will need one to begin.
I recommend getting python 3 from `anaconda <https://anaconda.org/>`_.
The Anaconda installer is straightforward to use, and you can let it set up your path, which makes the subsequent commands "just work."
However, Anaconda on Windows suggests not putting it on your path and instead using the Anaconda prompt.
The procedure is the same--just open an anaconda prompt window after installation then continue.
If you are on MacOS or Linux, you will want to restart your terminal after installing Anaconda so you get updated path specs.

Next, we need to install dependencies. GDAL is needed for accurate measurement of distance, and for converting coordinate systems.
I recommend getting it and segyio, used for interacting with the SEGY data format, using,

.. code-block:: bash

    conda install -c conda-forge gdal segyio

This step can be really slow, so don not worry if it is a bit painful.
At this point, I also recommend installing h5py--if you do not, it pip will do it on the next step. This can be done with

.. code-block:: bash

    conda install h5py

Now, you are ready to install impdar. You can get a version with

.. code-block:: bash

    pip install impdar

If you are not a super user, you may get an error related to permissions. This is fine, you just need to install for yourself only. Use 

.. code-block:: bash

    pip install --user impdar

You should now be all set to start using ImpDAR. Scroll down for documentation and links for examples.
    
Advanced
________

If you are not using Anaconda, you are on your own for installing dependencies. The challenges are generally GDAL and PyQt5, since these rely on libraries in other languages. For the most basic use cases, you can skip these, and go straight to installing ImpDAR with pip or through github.

To be sure that you have the newest version of ImpDAR as a lot of development is happening, you will want to use the development branch from GitHub. The pypi (pip) version is not updated as often to ensure a stable release. To get the devel version off git,

.. code-block:: bash

    git clone -b devel https://github.com/dlilien/ImpDAR.git
    cd impdar
    python setup.py install

This requires `git <https://git-scm.com/downloads>`_.

If you want to have the full suite of migration options, you will need to install `seisunix <https://github.com/JohnWStockwellJr/SeisUnix/wiki>`_.
The SeisUnix install is bit complicated, but there are instructions with it.
It should be possible to use SeisUnix on Windows with CygWin then interface with ImpDAR, but this is untested.
