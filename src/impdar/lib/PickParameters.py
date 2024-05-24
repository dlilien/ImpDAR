#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 David Lilien <dlilien90@gmail.com>
#
# Distributed under terms of the GNU GPL3.0 license.
"""Structure with input data for the picking algoriths."""


class PickParameters():
    """Some information used for determining for picks.

    This object contains several things that you need to know in order to
    pick a radar layer, like the frequency of layers you are looking for and
    the size window in which to search.

    Attributes
    ----------
    apickthresh: float
        Some kind of auto picking threshold (Unused: default 10)
    freq: float
        Frequency of the layer pick (default 4)
    dt: float
        Time between acquisitions
    plength: float
        The total packet to search for peaks
    FWW: float
        The width of the center portion which we are going to search
    scst: float
        The offset which we will search at
    pol: int
        Polarity of the picks
    apickflag: int
        I think this just kept track of whether StoDeep was autopicking
    addpicktype: str
        Some flag
    radardata: `RadarData`
        A link back up to the RadarData object with which this is affiliated
    """

    attrs = ['apickthresh',
             'freq',
             'dt',
             'plength',
             'FWW',
             'scst',
             'pol',
             'apickflag',
             'addpicktype']

    def __init__(self, radardata, pickparams_struct=None):
        if pickparams_struct is not None:
            for attr in self.attrs:
                setattr(self, attr, pickparams_struct[0][0][attr][0][0][0][0])
        else:
            self.freq = 4
            self.apickthresh = 10
            self.dt = radardata.dt
            self.pol = 1
            self.apickflag = 1
            self.addpicktype = 'zero'

        self.radardata = radardata
        self.freq_update(self.freq)

    def freq_update(self, freq):
        """Update the frequency at which we are looking.

        This is more complicated than just setting freq because other variables
        are a function of frequency and if not updated will break.

        Parameters
        ----------
        freq: float
            Target pick frequency.
        """
        self.freq = freq
        self.plength = 2 * int(round(1. / (
            self.freq * 1.0e6 * self.radardata.dt))) - 1
        if self.plength < 3:
            print('Warning: high freq compared to sampling rate. \
                  Forcing a minimum plength')
            self.plength = 3
        self.FWW = int(round(2. / 3. * (1. / (
            self.freq * 1.0e6 * self.radardata.dt))))
        if self.FWW % 2 == 0:
            self.FWW += 1
        self.scst = (self.plength - self.FWW) // 2

        # Guard against tiny datasets in this check...
        if self.plength > self.radardata.snum and self.radardata.snum >= 3:
            # print('Warning: Low freq compared to sampling rate. \
            #       Forcing a maximum plength')
            self.plength = self.radardata.snum
            self.FWW = self.radardata.snum // 2
            if self.FWW % 2 == 0:
                self.FWW += 1

    def to_struct(self):
        """Return attributes as a dictionary for saving.

        Guards against Nones so we can export to matlab

        Returns
        -------
        dict:
            The data for export with scipy.io.savemat
        """
        mat = {}
        for attr in self.attrs:
            if getattr(self, attr) is not None:
                mat[attr] = getattr(self, attr)
            else:
                mat[attr] = 0
        return mat
