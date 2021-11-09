#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 dlilien <dlilien@berens>
#
# Distributed under terms of the GNU GPL3.0 license.

"""This defines generic processing functions to ease calls from executables.

If interacting with the API, most processing steps should probably be called
by using methods on `RadarData` objects, so see
:doc:`that documentation<\./RadarData>` for most of your needs.
However, you may need to concatenate, which is defined separately because it
acts on multiple objects.

While the ``process`` and ``process_and_exit`` directives can be used, they
are generally not as useful as the direct calls.
"""
import os.path
import numpy as np

from .load import load
from .gpslib import interp as interpdeep
from .Picks import Picks

from copy import deepcopy


def process_and_exit(fn, cat=False, filetype='mat', o=None, **kwargs):
    """Perform one or more processing steps, save, and exit.

    Parameters
    ----------
    fn: list of strs
        The filename(s) to process.
    cat: bool, optional
        If True, concatenate files before processing rather than running
        through each individually.
    filetype: str, optional
        The type of input file. Default is .mat.
    o: str, optional
        An output path
    kwargs:
        These are the processing arguments for `process`
    """

    def _p_and_e(radar_data):
        processed = process(radar_data, **kwargs)
        if not processed and not cat:
            print('No processing steps performed. Not saving!')
        else:
            _save(radar_data, outpath=o, cat=cat)

    if cat:
        # first we do the quirky one
        # We need to load it all if catting
        radar_data = load(filetype, fn)
        radar_data = concat(radar_data)
        bn = os.path.splitext(fn[0])[0]
        if bn[-4:] == '_raw':
            bn = bn[:-4]
        radar_data[0].fn = bn + '_cat.mat'
        return _p_and_e(radar_data)
    else:
        # Otherwise, we can do things sequentially
        for fn_i in fn:
            radar_data = load(filetype, fn)
            return _p_and_e(radar_data)


def process(RadarDataList, interp=None, rev=False, vbp=None, hfilt=None,
            ahfilt=None, nmo=None, crop=None, hcrop=None, restack=None,
            denoise=None, migrate=None, **kwargs):
    """Perform one or more processing steps on a list of RadarData .

    Parameters
    ----------
    RadarDataList: list of strs
        The `~impdar.RadarData` objects to process
    rev: bool, optional
        Reverse the profile orientation. Default is False.
    vbp: 2-tuple, optional
        Vertical bandpass between (vbp1, vbp2) MHz.
        Default None (no filtering).
    hfilt: 2-tuple, optional
        Horizontal filter subtracting average trace between (hfilt1, hfilt2).
        Default is None (no hfilt).
    ahfilt: bool, optional
        Adaptively horizontally filter the data.
    denoise: bool, optional
        denoising filter (only wiener for now).
    migrate: string, optional
        Migrates the data.

    Returns
    -------
    processed: bool
        If True, we did something, if False we didn't
    """
    done_stuff = False

    # first some argument checking so we don't crash later
    if crop is not None:
        try:
            crop = (float(crop[0]), crop[1], crop[2])
        except ValueError:
            raise ValueError('First element of crop must be a float')
        except TypeError:
            raise TypeError('Crop must be subscriptible')
    if hcrop is not None:
        try:
            hcrop = (float(hcrop[0]), hcrop[1], hcrop[2])
        except ValueError:
            raise ValueError('First element of hcrop must be a float')
        except TypeError:
            raise TypeError('hcrop must be subscriptible')
        for dat in RadarDataList:
            dat.hcrop(*hcrop)
        done_stuff = True
    if denoise is not None:
        try:
            assert (type(denoise[0]) is int)
            assert (type(denoise[1]) is int)
        except (ValueError, TypeError, AssertionError, IndexError):
            raise ValueError('Denoise must be two integers giving vertical and horizontal window sizes')
    if vbp is not None:
        if not hasattr(vbp, '__iter__'):
            raise TypeError('vbp must be a tuple with first two elements \
                            [low] [high] MHz')
    if interp is not None:
        try:
            float(interp[0])
            interp[1]
        except (ValueError, TypeError, IndexError):
            raise ValueError('interp must be a target spacing (float) then a gps filename')


    if restack is not None:
        for dat in RadarDataList:
            if isinstance(restack, (list, tuple)):
                restack = int(restack[0])
            dat.restack(restack)
        done_stuff = True

    if rev:
        for dat in RadarDataList:
            dat.reverse()
        done_stuff = True

    if vbp is not None:
        for dat in RadarDataList:
            dat.vertical_band_pass(*vbp)
        done_stuff = True

    if hfilt is not None:
        for dat in RadarDataList:
            dat.hfilt(ftype='hfilt', bounds=hfilt)
        done_stuff = True

    if ahfilt:
        for dat in RadarDataList:
            dat.hfilt(ftype='adaptive',window_size=ahfilt)
        done_stuff = True

    if nmo is not None:
        if isinstance(nmo, (float, int)):
            print('One nmo value given. Assuming that this is the separation. \
                  Uice=1.6')
            nmo = (nmo, 1.6)
        for dat in RadarDataList:
            dat.nmo(*nmo)
        done_stuff = True

    if denoise is not None:
        for dat in RadarDataList:
            dat.denoise(*denoise)
        done_stuff = True

    if interp is not None:
        interpdeep(RadarDataList, float(interp[0]), interp[1])
        done_stuff = True

    # Crop after nmo so that we have nmo_depth available
    if crop is not None:
        for dat in RadarDataList:
            dat.crop(*crop)
        done_stuff = True

    if migrate is not None:
        for dat in RadarDataList:
            dat.migrate(mtype='stolt')
        done_stuff = True

    if not done_stuff:
        return False
    return True


def concat(radar_data):
    """Concatenate all radar data input.

    Parameters
    ----------
    radar_data: list of RadarData
        Objects to concatenate

    Returns
    -------
    RadarData
        A single, concatenated output.
    """
    # let's do some checks to make sure we are consistent here
    out = deepcopy(radar_data[0])

    for dat in radar_data[1:]:
        if out.snum != dat.snum:
            raise ValueError('Need the same number of samples in each file')
        if not np.allclose(out.travel_time, dat.travel_time):
            raise ValueError('Need matching travel time vectors')

    out.data = np.hstack([dat.data for dat in radar_data])
    tnums = np.hstack((np.array([0]),
                       np.cumsum([dat.tnum for dat in radar_data])))
    out.tnum = out.data.shape[1]
    out.trace_num = np.hstack([dat.trace_num + tnum for dat, tnum in zip(
                               radar_data, tnums)])
    if np.all([dat.dist is not None for dat in radar_data]):
        dists = np.hstack((np.array([0]),
                           np.cumsum([dat.dist[-1] for dat in radar_data])))
        out.dist = np.hstack([dat.dist + dist for dat, dist in zip(
            radar_data, dists)])
    for attr in ['pressure', 'trig', 'lat', 'long', 'x_coord', 'y_coord',
                 'elev', 'decday', 'trace_int']:
        if np.all([getattr(dat, attr) is not None for dat in radar_data]):
            setattr(out, attr, np.hstack([getattr(dat, attr)
                                          for dat in radar_data]))

    # Picks are the most challenging part
    all_picks = []
    for dat in radar_data:
        if dat.picks is not None and dat.picks.picknums is not None and dat.picks.picknums != 0:
            all_picks.extend(dat.picks.picknums)
    all_picks = np.unique(all_picks).tolist()
    out.picks = Picks(out)
    if len(all_picks) > 0:
        out.picks.picknums = all_picks
        out.picks.lasttrace.tnum = [out.tnum for i in all_picks]
        out.picks.lasttrace.snum = [0 for i in all_picks]
        pick_attrs = ['samp1', 'samp2', 'samp3', 'power', 'time']
        for attr in pick_attrs:
            setattr(out.picks, attr, np.zeros((len(all_picks), out.tnum)) * np.NaN)
        start_ind = 0
        for dat in radar_data:
            if ((not hasattr(dat, 'picks')) or (not hasattr(dat.picks, 'picknums')) or (
                    not hasattr(dat.picks.picknums, '__len__')) or (len(dat.picks.picknums) == 0)):
                start_ind += dat.tnum
                continue
            for attr in pick_attrs:
                if hasattr(dat.picks, attr):
                    in_dat = getattr(dat.picks, attr)
                    if in_dat is not None:
                        out_dat = getattr(out.picks, attr)
                        for pick in dat.picks.picknums:
                            out_dat[all_picks.index(pick), start_ind:start_ind + dat.tnum] = in_dat[
                                dat.picks.picknums.index(pick), :]
                        setattr(out.picks, attr, out_dat)
            start_ind += dat.tnum

    print('Objects concatenated')
    return [out]


def _save(rd_list, outpath=True, cat=False):
    if outpath is not None:
        if len(rd_list) > 1:
            for rd in rd_list:
                bn = os.path.split(os.path.splitext(rd.fn)[0])[1]
                if bn[-4:] == '_raw':
                    bn = bn[:-4]
                out_fn = os.path.join(outpath, bn + '_proc.mat')
                rd.save(out_fn)
        else:
            out_fn = outpath
            rd_list[0].save(out_fn)
    else:
        for rd in rd_list:
            bn = os.path.splitext(rd.fn)[0]
            if bn[-4:] == '_raw':
                bn = bn[:-4]
            if cat:
                out_fn = bn + '.mat'
            else:
                out_fn = bn + '_proc.mat'
            rd.save(out_fn)
