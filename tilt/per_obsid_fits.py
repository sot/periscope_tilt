import os
import logging
import numpy as np

import cPickle
import re

# Matplotlib setup                                                                                                     
# Use Agg backend for command-line (non-interactive) operation                                                         
import matplotlib
if __name__ == '__main__':
    matplotlib.use('Agg')
import matplotlib.pyplot as plt


import sherpa.ui as ui

from Ska.Numpy import interpolate


logger = logging.getLogger("fit")

DATADIR = 'auto'
USER_PARS = ['axial', 'diam']


def tilt_model(tilt_data, evt_times, user_pars=None):
    if user_pars is None:
        user_pars = USER_PARS
    def model(pars, x):
        tilt_times = tilt_data['time']
        tilt_axial = tilt_data['tilt_axial']
        tilt_diam = tilt_data['tilt_diam']
        ax = interpolate(tilt_axial, tilt_times , evt_times)
        di = interpolate(tilt_diam, tilt_times, evt_times)
        if len(user_pars) == 4:
            m = (pars[0] * 3600 * (ax + pars[2])
                 + pars[1] * 3600 * (di + pars[3]))
        elif len(user_pars) == 2:
            m = (pars[0] * 3600 * ax + pars[1] * 3600 * di)
        else:
            raise ValueError("have {} user pars".format(len(user_pars)))
        return m - np.mean(m)
    return model


def binned_mean(data, evtstime, tbin=5000.):
    ts = []
    ds = []
    dsminus = []
    dsplus = []
    for win_start in range(0,
                           int(evtstime[-1]-evtstime[0])-int(tbin),
                           int(tbin)):
        tmask = ((evtstime-evtstime[0] >= win_start)
                 & (evtstime-evtstime[0] < win_start + int(tbin)))
        range_data = data[ tmask ]
        range_time = evtstime[ tmask ]
        if np.std(range_data) > 0:
            ds.append(np.mean(range_data))
            dsminus.append(np.std(range_data)/np.sqrt(len(range_data)))
            dsplus.append(np.std(range_data)/np.sqrt(len(range_data)))
            t = np.mean(range_time)
            ts.append(t)
    return [ ts, ds, dsminus, dsplus ]



def plot_fits(ids, posdir=None, outdir='fit'):

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    for id in ids:
        obsid = id['obsid']
        data_id = id['data_id']
        ax = id['ax']
        pf = open(os.path.join(posdir, 'pos.pkl'), 'r')
        pos = cPickle.load(pf)
        pf.close()

        pos_data = pos[ax]

        fig = plt.figure(figsize=(6,4))
        axplot = ui.get_model_plot(data_id)
        plt.plot( axplot.x/1000., axplot.y, 'b.')
        [ts, ds, dsminus, dsplus] = binned_mean(
            pos_data-np.mean(pos_data),
            pos['time']-pos['time'][0] )
        if not len(ts):
            continue
        plt.errorbar((ts-ts[0])/1000., ds,
                     yerr=dsminus, color='k', marker='.', linestyle='')
        plt.xlabel('Time (ksec)')
        plt.ylabel('%s' % ax)
        plt.subplots_adjust(left=.15)
        plt.title('Model(Blue), BinMean Pos(Black) (obs=%s)' % (obsid))
        plt.ylim(-.3,.3)
        plt.grid()
        plt.savefig('%s/%s_%s_sfit.png' % (outdir, obsid, ax))
        plt.close()



def run_fits(obsids, ax, user_pars=None,
             fixed_pars=None, guess_pars=None, label='model',
             per_obs_dir='per_obs_nfits',
             outdir=None, redo=False):

    if len(obsids) == 0:
        print "No obsids, nothing to fit"
        return None
    if user_pars is None:
        user_pars = USER_PARS

    if not os.path.exists(per_obs_dir):
        os.makedirs(per_obs_dir)

    obsfits = []
    for obsid in obsids:

        outdir = os.path.join(per_obs_dir, 'obs{:05d}'.format(obsid))
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        model_file = os.path.join(outdir, '{}.pkl'.format(label))
        if os.path.exists(model_file) and not redo:
            #logger.warn('Using previous fit found in %s' % model_file)
            print model_file
            mod_pick = open(model_file, 'r')
            modelfit = cPickle.load( mod_pick )
            mod_pick.close()
            obsfits.append(modelfit)
            continue

        modelfit = {'label': obsid}

        ui.clean()
        data_id = 0
        obsdir = "%s/obs%05d" % (DATADIR, obsid)
        tf = open(os.path.join(obsdir,'tilt.pkl'), 'r')
        tilt = cPickle.load(tf)
        tf.close()
        pf = open(os.path.join(obsdir, 'pos.pkl'), 'r')
        pos = cPickle.load(pf)
        pf.close()

        pos_data = pos[ax]
        point_error = 5
        pos_data_mean = np.mean(pos_data)
        ui.set_method('simplex')

        # Fit a line to get more reasonable errors
        init_staterror = np.zeros(len(pos_data))+point_error
        ui.load_arrays(data_id,
                       pos['time']-pos['time'][0],
                       pos_data-np.mean(pos_data),
                       init_staterror)
        ui.polynom1d.ypoly
        ui.set_model(data_id, 'ypoly')
        ui.thaw(ypoly.c0, ypoly.c1)
        ui.fit(data_id)
        fit = ui.get_fit_results()
        calc_staterror = init_staterror * np.sqrt(fit.rstat)
        ui.set_staterror(data_id, calc_staterror)
        # Confirm those errors
        ui.fit(data_id)
        fit = ui.get_fit_results()
        if ( abs(fit.rstat-1) > .2):
            raise ValueError('Reduced statistic not close to 1 for error calc')

        # Load up data to do the real model fit
        fit_times = pos['time']
        tm_func = tilt_model(tilt,
                             fit_times,
                             user_pars=user_pars)

        ui.get_data(data_id).name = str(obsid)
        ui.load_user_model(tm_func, 'tiltm%d' % data_id)
        ui.add_user_pars('tiltm%d' % data_id, user_pars)
        ui.set_method('simplex')
        ui.set_model(data_id, 'tiltm%d' % (data_id))
        ui.set_par('tiltm%d.diam' % data_id, 0)

        if fixed_pars is not None and ax in fixed_pars:
            for par in fixed_pars[ax]:
                ui.set_par('tiltm{}.{}'.format(0, par), fixed_pars[ax][par])
                ui.freeze('tiltm{}.{}'.format(0, par))

        if guess_pars is not None and ax in guess_pars:
            for par in guess_pars[ax]:
                ui.set_par('tiltm{}.{}'.format(0, par), guess_pars[ax][par])

        ui.show_all()
        # Fit the tilt model
        ui.fit(data_id)
        fitres = ui.get_fit_results()
        ui.confidence(data_id)
        myconf = ui.get_confidence_results()

#        save_fits(ax=ax, fit=fitres, conf=myconf, outdir=outdir)
#        plot_fits(ids,outdir=os.path.join(outdir,'fit_plots'))

        axmod = dict(fit=fitres, conf=myconf)
        for idx, modpar in enumerate(myconf.parnames):
            par = modpar.lstrip('tiltm0.')
            axmod[par] = ui.get_par('tiltm0.%s' % par).val
            axmod["{}_parmax".format(par)] = myconf.parmaxes[idx]
            axmod["{}_parmin".format(par)] = myconf.parmins[idx]
        modelfit[ax] = axmod

        mod_pick = open(model_file, 'w')
        cPickle.dump( modelfit, mod_pick)
        mod_pick.close()

        obsfits.append(modelfit)

        plot_fits([dict(obsid=obsid, data_id=data_id, ax=ax)],
                  posdir=obsdir,
                  outdir=outdir)


    return obsfits

