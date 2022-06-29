Apres Radar (apdar) examples
===================

``apdar`` is the ApRES (Autonomous phase-sensitive Radio Echo Sounder) module for ImpDAR. Most of the source code originally came from Matlab scripts written at BAS, and credit is given to those scripts where appropriate.

apdar load
-------

With ``apdar load``, you can import an ApRES data file from one of several formats (.mat, .hdf5, .DAT, .dat). The load scripts will determine the data type for you and assign appropriate fields to an apres data object. The raw data do not mean much, but you can make sure that some data are written to the file at least by plotting, ``apdar plot apres_raw.mat``,

.. image:: apres_raw.png

apdar processing
-------

``apdar range 4000 apres_raw.mat``

This does a range conversion (pulse compression) using the known transmit chirp which is written into the file header. Now we have range and can plot amplitude or power against it, more like a traditional radar A-scope image, ``apdar plot apres_range.mat``,

.. image:: apres_range.png

This is only plotting the first chirp whereas the file will have many chirps in each burst and typically many bursts as well. To stack the chirps use,

``apdar stack 0 apres_range.mat``

 Prescribing '0' stacks *all* the chirps into one, including across bursts. You can also give a specific number of chirps to stack.

.. image:: apres_range_stacked.png

Calculate the complex uncertainty for each sample using,

``apdar uncertainty 3000 apres_range_stacked.mat``

This looks at samples below a given range (3000 given here) to determine a 'median noise phasor' which is used to calculate the uncertainty for each sample by comparing them to the noise phasor.

.. image:: apres_range_stacked_uncertainty.png


apdar diffload
-------

With ``apdar diffload apres_1.mat apres_2.mat``, you can load two ApRES acquisitions side-by-side in order to difference them for a vertical velocity calculation. This only takes .mat or .hdf5 files as input, and those files should be stacked and range converted beforehand.

.. image:: diff_raw.png

apdar difference processing
-------

``apdar pdiff apresdiff_raw.mat``

This command calculates the phase coherence between the two acquisitions.

.. image:: diff_coherence.png

``apdar diffproc apresdiff_raw.mat``

This command walks through the entire processing flow for differencing two ApRES acquisitions.

.. image:: diff_velocity.png
