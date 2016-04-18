import os
import sys
import math
import logging
import numpy as np
import time
from glob import glob
import cPickle
import re

# Matplotlib setup
# Use Agg backend for command-line (non-interactive) operation
import matplotlib
if __name__ == '__main__':
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

import sherpa.ui as ui
import scipy.special

from Chandra.Time import DateTime
from Ska.Shell import bash
import Ska.DBI
from Ska.Numpy import interpolate


logger = logging.getLogger("fit")
DATADIR = 'auto'


def tilt_model(tilt_data, evt_times):
    tilt_times = tilt_data['time']
    tilt_axial = tilt_data['tilt_axial']
    tilt_diam = tilt_data['tilt_diam']
    ax = interpolate(tilt_axial, tilt_times , evt_times)
    di = interpolate(tilt_diam, tilt_times, evt_times)
    bulk = interpolate(tilt_data['OHRTHR42'].vals, tilt_data['OHRTHR42'].times, evt_times)
    def model(pars, x):
        m = (pars[0] * 3600 * ax  + pars[1] * 3600 * di)
        return m - np.mean(m)
    return model



def get_options():
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_defaults()
    parser.add_option("--outdir",
                      default="rate_limit",
                      help="Output directory")
    parser.add_option("--fit_all_pts",
                      default=True,
                      action='store_true',
                      help="Fit every telemetry point")
    parser.add_option("--obsid",
                      action='append',
                      type=int)
    opt, args = parser.parse_args()
    return opt, args


def setup_fit(obsids, tilt_model, ax='yag', user_pars=['axial', 'diam', ], opt=None):

    ids = []

    for data_id, obsid in enumerate(obsids):
        obsdir = "%s/obs%05d" % (DATADIR, obsid)
        tf = open(os.path.join(obsdir,'tilt.pkl'), 'r')
        tilt = cPickle.load(tf)
        tf.close()
        pf = open(os.path.join(obsdir, 'pos.pkl'), 'r')
        pos = cPickle.load(pf)
        pf.close()

        point_error = 5
        ui.set_method('simplex')
        fit_times = None

        if opt.fit_all_pts:
            all_time = pos['time']
            # limit to a count rate of 0.2 counts/sec
            count_rate = 0.2
            obs_dur = all_time[-1] - all_time[0]
            counts = obs_dur * count_rate
            # sample with replacement
            index = np.random.choice(np.arange(len(pos[ax])), counts, replace=True)
            pos_data = pos[ax][index]
            time = pos['time'][index]

            init_staterror = np.zeros(counts)+point_error
            ui.load_arrays(data_id,
                           time-pos['time'][0],
                           pos_data - np.mean(pos[ax]),
                           init_staterror)
            ui.polynom1d.ypoly
            ui.set_model(data_id, 'ypoly')
            ui.thaw(ypoly.c0, ypoly.c1)

            ui.fit(data_id)
            fit = ui.get_fit_results()
            calc_staterror = init_staterror * np.sqrt(fit.rstat)
            ui.set_staterror(data_id, calc_staterror)
            ui.fit(data_id)
            fit = ui.get_fit_results()
            if ( abs(fit.rstat-1) > .2):
                raise ValueError('Reduced statistic not close to 1 for error calc')
            fit_times = time
        else:
            [ ts, ds, dsminus, dsplus] = binned_mean(
                pos_data-np.mean(pos_data),
                pos['time'])
            fit_times = ts
            ui.load_arrays(data_id,
                           ts-ts[0],
                           ds-np.mean(ds),
                           dsminus)

        tm_func = tilt_model(tilt,
                             fit_times)

        ui.get_data(data_id).name = str(obsid)
        model_label = 'tiltm{}'.format(data_id)
        ui.load_user_model(tm_func, model_label)
        ui.add_user_pars(model_label, user_pars)
        ui.set_method('simplex')
        ui.set_model(data_id, model_label)
        ids.append(dict(obsid=obsid,
                        data_id=data_id,
                        ax=ax))
    return ids

def save_fits(ax, fit, conf=None, outdir='fit_nodiam'):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    fit_file = open('%s/%s_fit.txt' % (outdir, ax), 'w')
    fit_file.write(str(fit) + "\n" + str(conf) + "\n")
    fit_file.close()



def fit(obsids, opt, tilt_model, user_pars=['axial','diam'],
        myguess=None, redo=False):

    outdir = opt.outdir

    model_file = os.path.join(outdir, 'both_float_model.pkl')
    if os.path.exists(model_file) and not redo:
        logger.warn('Using previous fit found in %s' % model_file)
        mod_pick = open(model_file, 'r')
        modelfit = cPickle.load( mod_pick )
        mod_pick.close()
        return modelfit

    idlist = []
    modelfit = {'label': obsids[0]}

    for ax in ['yag', 'zag']:
        ui.clean()
                  
        ids = setup_fit(obsids, tilt_model,
                        ax=ax, user_pars=user_pars, opt=opt)

        for id in ids:
            if id['data_id'] != 0:
                for par in user_pars:
                    ui.link('tiltm%d.%s' % (id['data_id'], par),
                            'tiltm%d.%s' % (0, par))
            idlist.append(id['data_id'])

        if myguess:
            for par in myguess[ax].keys():
                ui.set_par('tiltm0.%s' % par, myguess[ax][par])

        idlist.sort()
        ui.fit(*idlist)
        ui.freeze('tiltm0.diam')
        fitres = ui.get_fit_results()
        ui.thaw('tiltm0.diam')
        ui.fit(*idlist)
        fitres = ui.get_fit_results()
        ui.confidence(*idlist)
        myconf = ui.get_confidence_results()

        save_fits(ax=ax, fit=fitres, conf=myconf, outdir=outdir)
        #plot_fits(ids,outdir=os.path.join(outdir,'fit_plots'))

        axmod = dict(fit=fitres,conf=myconf)
        for idx, modpar in enumerate(myconf.parnames):
            par = modpar.lstrip('tiltm0.')
            axmod[par] = ui.get_par('tiltm0.%s' % par).val
            axmod["{}_parmax".format(par)] = myconf.parmaxes[idx]
            axmod["{}_parmin".format(par)] = myconf.parmins[idx]
        modelfit[ax] = axmod


    mod_pick = open(os.path.join(outdir, 'both_float_model.pkl'), 'w')
    cPickle.dump( modelfit, mod_pick ) 
    mod_pick.close()
    return modelfit




if __name__ == '__main__':

    opt, args = get_options()

    from astropy.table import Table
    srcs = Table.read('src_table.dat', format='ascii')
    obsids = srcs[(srcs['NET_COUNTS'] > 2500) & (srcs['point_source'] == 'True')]['obsid']
    #obsids = [15803, 17128, 15192, 17174,  17561]
    #obsids = [10764]
    #               2016    2015   2014   2013   202    2013   2012   2011   2010,  2009   2008  2006
    #if opt.obsid is not None:
    #    obsids = opt.obsid

    #obsids = [17105, 17128, 14567, 14457, 14567, 13305, 10764]

    user_pars = ['axial', 'diam']
    modelfit = fit(obsids, opt, tilt_model, user_pars)
#        resid_obsids = [17105, 17128, 16221, 15741, 14567, 14457, 14567, 13305, 12809, 10062, 10764, 7302]
#
#    resids = get_resids(obsids,
#                        tilt_model, modelfit, user_pars,
#                        opt.outdir)
#    resid_limits = dict(abs=.15,
#                        rel=5)
#    plot_resids(resids,
#                os.path.join(opt.outdir, 'resid_plots'),
#                resid_limits)
#    index(
#        modelfit,
#        resids,
#        opt.outdir,
#        'resid_plots',
#        resid_limits
#        )
#
