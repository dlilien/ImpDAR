Interpretation
==============

Interpretation in this context primarily means picking layers (either isochrones or the bed). In the future, this functionality may be expanded to make picking other things, e.g. some discontinuity, easier.

Functions used for picking
--------------------------

.. automodule:: impdar.lib.picklib
    :members:


Classes used by interpreter
---------------------------

These classes are broken down to match the structure of StODeep, so we store information about how the picks get made, and the picks themselves, using different objects.

If you have done some interpretation, you will likely want to subsequently interact with the `Picks` object. Often, this can be done without accessing the API by converting the picks/geospatial data to another form, e.g. with `impdar convert shp fn_picked.mat`. You can also make plots with the picks on top of the radar data, or with the return power in geospatial coordinates, using `impplot rg fn_picked.mat` or `impplot power fn_picked.mat layer_num`. For further operations, you will probably want to access the `Picks` object described next. For example, using the picks object you could do something like

.. code-block:: python

    import numy as np
    import matplotlib.pyplot as plt
    from impdar.lib import RadarData
    rd = RadarData('[PICKED_DATA_FN.mat]')

    # make a basic plot of the radargram
    fig, ax = plt.subplots()
    im, _, _, _, _ = plot.plot_radargram(rd, fig=fig, ax=ax, xdat='dist', ydat='depth', return_plotinfo=True)

    # calculate the return power
    c = 10. * np.log10(rd.picks.power[0, :])
    c -= np.nanmax(c)

    # plot the return power on the layer, being careful of NaNs
    mask = ~np.isnan(rd.picks.samp1[0, :])
    cm = ax.scatter(rd.dist[mask.flatten()],
                    rd.nmo_depth[rd.picks.samp1[0, :].astype(int)[mask]],
                    c=c.flatten()[mask.flatten()],
                    s=1)


.. automodule:: impdar.lib.Picks
    :members:

.. automodule:: impdar.lib.PickParameters
    :members:
