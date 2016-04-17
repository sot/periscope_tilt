import os
import numpy as np
from astropy.table import Table, vstack
import cPickle
from Chandra.Time import DateTime
from Ska.Numpy import interpolate
from itertools import izip


def binned_mean(data, evtstime, tbin=1000.):
    """
    bin data in time bins of tbin size using evtstime as time data
    :param data: data to bin as array or numpy array
    :param evtstime: times associated with data
    :param tbin: bin time in units that match evtstime (generally seconds)
    :returns: [bin times,
                data mean for time bin,
                data std of bin]
    """
    ts = []
    ds = []
    ds_std = []
    for win_start in range(0,
                           int(evtstime[-1]-evtstime[0])-int(tbin),
                           int(tbin)):
        tmask = ((evtstime-evtstime[0] >= win_start)
                 & (evtstime-evtstime[0] < win_start + int(tbin)))
        range_data = data[ tmask ]
        range_time = evtstime[ tmask ]
        if np.std(range_data) > 0:
            ds.append(np.mean(range_data))
            ds_std.append(np.std(range_data)/np.sqrt(len(range_data)))
            #t = np.mean(range_time)
            t = (np.max(range_time) + np.min(range_time))/2.
            ts.append(t)
    return [ts, ds, ds_std]



def per_obs_binned_data(obs):
    """
    Fetch X-ray position and gradient telemetry for a given observation,
    bin it by time, and return a brief table (astropy Table) of useful values where each row
    corresponds to a bin time slice.  This table includes:

    time (bin mean time),
    yag, y_std, yag_mean (for whole observation not the bin),
    zag, z_std, zag_mean (for whole observation not the bin),
    oobagrd3, oobagrd3_mean (for observation),
    oobagrd6, oobagrd6_mean (for observation),
    ohrthr42, ohrthr42_mean (for observation),
    obsid

    The "mean" and obsid cols are fixed per observation, and are just conveniences to work with the
    time binned data
            """
    obsdir = os.path.join(DATADIR, 'obs{:05d}'.format(obs))
    pos_data = cPickle.load(open(os.path.join(obsdir, 'pos.pkl')))
    tilt = cPickle.load(open(os.path.join(obsdir, 'tilt.pkl')))
    ts, yag, yag_std = binned_mean(pos_data['yag'], pos_data['time'], BIN)
    ts, zag, zag_std = binned_mean(pos_data['zag'], pos_data['time'], BIN)
    obs_chunk = {'time': ts, 'yag': yag, 'yag_std': yag_std, 'zag': zag, 'zag_std': zag_std}
    obs_chunk['yag_mean'] = np.repeat(np.mean(pos_data['yag']), len(ts))
    obs_chunk['zag_mean'] = np.repeat(np.mean(pos_data['zag']), len(ts))
    for msid in ['oobagrd3', 'oobagrd6', 'ohrthr42']:
        obs_chunk[msid] = interpolate(tilt[msid.upper()].vals, tilt[msid.upper()].times, ts)
        ok_time = ((tilt[msid.upper()].times > pos_data['time'][0]) &
                   (tilt[msid.upper()].times < pos_data['time'][-1]))
        obs_chunk['{}_mean'.format(msid)] = np.repeat(np.mean(tilt[msid.upper()].vals[ok_time]), len(ts))
    obs_chunk['obsid'] = np.repeat(obs, len(ts)).tolist()
    return Table(obs_chunk)


def get_all_bins(resid_obsids):
    """
    Get binned data for all requested obsids and make one big astropy table.  
    See per_obs_binned_data for a brief description of table columns
    """
    bin_data = None
    for obsid in resid_obsids:
        odata = per_obs_binned_data(obsid)
        if bin_data is None:
            bin_data = odata
        else:
            bin_data = vstack([bin_data, odata])
    return bin_data


def get_obs_col_p2p(obsid, col, table):
    """
    Given an astropy table from get_all_bins, a specific obsid, and a requested column, return the
    min and max values of the column during the obsid in a dictionary
    """
    obs_data = table[table['obsid'] == obsid]
    return {'p2p_{}'.format(col): np.max(obs_data[col]) - np.min(obs_data[col]),
            'obsid': obsid, 'time': obs_data[0]['time']}



def get_obs_col_std(obsid, col, table):
    """
    Given an astropy table from get_all_bins, a specific obsid, and a requested column, return the
    std of the column during the obsid in a dictionary
    """
    obs_data = table[table['obsid'] == obsid]
    return {'std_{}'.format(col): np.std(obs_data[col]), 'obsid': obsid, 'time': obs_data[0]['time']}


