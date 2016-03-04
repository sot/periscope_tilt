"""
Walks through all science observations, retrieves the src2 file for each, finds the source with the most
counts from that file, and adds the source to the 'bright_source' table in the dt.db3 sqlite3 database
"""
import os
import numpy as np
from glob import glob
import json
import tempfile
from datetime import datetime
import cPickle

from astropy.table import Table
from astropy.io import fits
from Ska.Numpy import smooth, interpolate
from Ska.Shell import bash, tcsh_shell, getenv
from Ska.DBI import DBI
from Ska.engarchive import fetch
import Quaternion
import Ska.quatutil




REDO = False

sqlaca = DBI(dbi='sybase', user='aca_read')



#XRAY_DATA = '/data/aca/archive/xray_for_periscope'
projdir = '/proj/sot/ska/analysis/periscope_tilt_2016'
XRAY_DATA = os.path.join(projdir, 'auto')
ciao_env = getenv("source /soft/ciao/bin/ciao.csh", shell='tcsh')


def get_on_axis_bright(srctable, x_center, y_center, limit=180):
    """
    Given source table from celldetect and the X/Y center for the observation
    return a reduced source list with just the brightest source within radius limit

    :param srctable: table of sources containing at least X, Y, and NET_COUNTS for each source
    :param x_center: rough aimpoint X center from evt2 file
    :param y_center: rough aimpoint Y center from evt2 file
    :param limit: radial limit in pixels
    :returns: reduced source table containing just brightest source within radius limit
    """
    dist = np.sqrt(((srctable['X'] - x_center)**2
                    + (srctable['Y'] - y_center)**2))
    # limit to sources within ~90 arcsecs
    if not np.any(dist < limit):
        return None
    srctable = srctable[dist < limit]
    return srctable[np.argmax(srctable['NET_COUNTS'])]


def get_source(obsid, obi, obsdir):
    """
     Get bright source for an observation from previously fetched data

    :param obs: observations table entry for obsid
    :param obsdir: local directory for observation
    :returns: astropy table with one row describing bright source
    """
    srcfiles = glob(os.path.join(obsdir, '*src2*'))
    if len(srcfiles) == 0:
        return None
    src_json = '{}/evtinfo.json'.format(obsdir)
    proc_info = json.load(open(src_json))
    if int(proc_info['OBS_ID']) != obsid:
        raise ValueError("Unexpected obsid mismatch")
    maxsrc = None
    srcfile = srcfiles[0]
    srctable = Table.read(srcfile)
    if len(srctable):
        maxsrc = get_on_axis_bright(srctable, proc_info['x_center'], proc_info['y_center'])
    if not maxsrc:
        return None
    maxsrc = Table(maxsrc)[['X', 'Y', 'RA', 'DEC', 'NET_COUNTS', 'SNR', 'DETSIZE']]
    maxsrc['obsid'] = obsid
    maxsrc['obi'] = obi
    maxsrc['ascdsver'] = proc_info['ASCDSVER']
    maxsrc['caldbver'] = proc_info['CALDBVER']
    maxsrc['revision'] = proc_info['REVISION']
    return maxsrc


def find_obsid_src(obsid, obs):
    """
    For an obsid, fetch evt2 data as needed, run celldetect, and return information about
    a bright source if available.

    :param obsid: obsid
    :returns:  astropy table with row describing bright source

    """
    if obs['readmode'] == 'CONTINUOUS':
        return None

    xray_data = XRAY_DATA
    if not os.path.exists(xray_data):
        os.makedirs(xray_data)
    obsdir = os.path.join(xray_data, 'obs%05d' % obsid)
    if not os.path.exists(obsdir):
        os.makedirs(obsdir)

    no_source = os.path.join(obsdir, "has_no_source")
    if os.path.exists(no_source):
        return None

    srcfiles = glob(os.path.join(obsdir, '*src2*'))
    src_json = '{}/evtinfo.json'.format(obsdir)
    if len(srcfiles) and os.path.exists(src_json) and not REDO :
        src = get_source(obs['obsid'], obs['obi'], obsdir)
        return src

    # Fetch X-ray events
    tempdir = tempfile.mkdtemp(dir='/export/jeanconn/tempdir/')
    det = 'acis'
    if (obs['instrume'] == 'HRC'):
        det = 'hrc'
    bash('echo "cd %s\n obsid=%d\n get %s2{evt2}\n" | arc5gl -stdin' % (tempdir, obs['obsid'], det))
    event_files = glob("%s/*evt2*gz" % tempdir)
    if not len(event_files):
        return None
    if len(event_files) > 1:
        return None
    bash('gunzip %s' % event_files[0])
    event_files = glob("%s/*evt2*" % tempdir)
    # Read processing keywords from the evt2 file
    f = fits.open(event_files[0])
    proc_info = {'ASCDSVER': f[1].header['ASCDSVER'],
                 'CALDBVER': f[1].header['CALDBVER'],
                 'REVISION': f[1].header['REVISION'],
                 'ASOLFILE': f[1].header['ASOLFILE'],
                 'OBS_ID': f[1].header['OBS_ID']}
    # Get the X/Y for this detector from the header
    for i in range(1, 30):
        if 'TTYPE{}'.format(i) in f[1].header:
            if f[1].header['TTYPE{}'.format(i)] == 'x':
                proc_info['x_center'] = f[1].header['TCRPX{}'.format(i)]
                proc_info['y_center'] = f[1].header['TCRPX{}'.format(i + 1)]
            if f[1].header['TTYPE{}'.format(i)] == 'detx':
                proc_info['detx_center'] = f[1].header['TCRPX{}'.format(i)]
                proc_info['dety_center'] = f[1].header['TCRPX{}'.format(i + 1)]
    if ((proc_info['detx_center'] != proc_info['x_center'])
        or (proc_info['dety_center'] != proc_info['y_center'])):
        raise ValueError
    f.close()
    # Save processing keywords out to a file
    ph = open(src_json, 'w')
    ph.write(json.dumps(proc_info, indent=4, sort_keys=True))
    ph.close()
    # Setup to run celldetect
    pfiles = ';'.join([tempdir, ciao_env['PFILES'].split(';')[-1]])
    outlines1, stat = tcsh_shell("env PFILES={} punlearn celldetect".format(pfiles), env=ciao_env)
    cmd = ('{} celldetect eenergy=0.9 infile="{}" outfile="{}/{}_src2.fits"'.format(
            'env LD_LIBRARY_PATH="" PFILES="{}"'.format(pfiles), event_files[0], obsdir, obs['obsid']))
    # Run celldetect
    outlines2, stat = tcsh_shell(cmd, env=ciao_env)
    srcfiles = glob(os.path.join(obsdir, "*src2*"))
    if len(srcfiles):
        src = get_source(obs['obsid'], obs['obi'], obsdir)
        if src is not None:
            os.unlink(event_files[0])
        else:
            bash("touch {}".format(no_source))
            print "no sources in {}".format(event_files[0])
        return src
    else:
        raise ValueError("No src2 file was made")



def filter_bad_telem(msid, method='nearest'):
    # use the bad quality field to select
    # and replace bad data in place using the given method
    ok = msid.bads == False
    bad = msid.bads == True
    fix_vals = interpolate(msid.vals[ok],
                           msid.times[ok],
                           msid.times[bad],
                           method=method)
    msid.vals[bad] = fix_vals
    return msid


def extract_point(obs_info, src, obsdir, point):
    det = 'acis'
    if obs_info['instrume'] == 'HRC':
        det = 'hrc'
    tempdir = tempfile.mkdtemp()
    bash('echo "cd %s\n obsid=%d\n get %s2{evt2}\n" | arc5gl -stdin' %
         (tempdir, src['obsid'], det))
    reg = os.path.join(obsdir, 'center.reg')
    c = open(reg, 'w')
    regstring = "circle(%f, %f, 6)" % (src[0]['X'], src[0]['Y'])
    c.write("%s\n" % regstring)
    c.close()
    evt2 = glob('%s/*evt2.fits*' % tempdir)[0]
    dmstring = '[cols time,ra,dec,x,y]'
    if det == 'acis':
        dmstring = dmstring +  '[energy=300:7000]'
    #print("/proj/sot/ska/bin/doapp -ciao dmcopy %s'[(x,y)=%s]%s' %s" %
    #      (evt2, regstring, dmstring, point))
    status = bash("/proj/sot/ska/bin/doapp -ciao dmcopy %s'[(x,y)=%s]%s' %s" %
                  (evt2, regstring, dmstring, point))
    if not status:
        os.unlink(evt2)


def get_xray_data(obsids=None):

    if obsids is None:
        obsid_query = "select obsid from obs_periscope_tilt where max_oobagrd3 - min_oobagrd3 > .09 order by (max_oobagrd3 - min_oobagrd3) desc"
        obsids = sqlaca.fetchall(obsid_query)['obsid']
        print "Fetching obsids using query '{}'".format(obsid_query)
        print "Looking for sources in these {} observations".format(len(obsids))
    else:
        print "Fetching specified obsids"
        print obsids

    for obsid in obsids:

        obsdir = "%s/auto/obs%05d" % (projdir, obsid)
        if not os.path.exists(obsdir):
            os.makedirs(obsdir)

        src_file = os.path.join(obsdir, 'picked_src.dat')
        obs_info = sqlaca.fetchone("select * from observations where obsid = %d" % obsid)
        if not os.path.exists(src_file) or REDO:
            src = find_obsid_src(obsid, obs_info)
            if src is None:
                continue
            else:
                src.write(src_file,
                          format='ascii.tab')
        else:
            src = Table.read(src_file, format='ascii.tab')


        # cut-out region with a point source for this obsid
        point = '%s/point_source.fits' % obsdir
        #print point
        if not os.path.exists(point) or REDO:
            extract_point(obs_info, src, obsdir, point)
        obsid_src = point

        # periscope tilt telemetry
        if not os.path.exists(os.path.join(obsdir, 'tilt.pkl')) or REDO:
            obs = obs_info
            msids=[ 'OOBAGRD3', 'OOBAGRD6', 'OHRTHR42', 'OHRTHR43', 'OOBTHR39', 'OHRTHR24', '4RT702T', 'OHRTHR24', 'AACBPPT', 'AACH1T', 'AACCCDPT', 'AACBPRT']
            telem = fetch.MSIDset(msids, obs['tstart'], obs['tstop'])
            telemtime = telem['OOBAGRD3'].times
            tilt = dict(telem)
            tilt.update(dict(
                time = telemtime,
                tilt_axial = telem['OOBAGRD3'].vals,
                tilt_diam = telem['OOBAGRD6'].vals))
            tilt_pick = open(os.path.join(obsdir, 'tilt.pkl'), 'w')
            cPickle.dump( tilt, tilt_pick)
            tilt_pick.close()
            #print os.path.join(obsdir, 'tilt.pkl')


        # position data
        if not os.path.exists(os.path.join(obsdir, 'released_pos.pkl')) or REDO:
            obs = obs_info
            evts = Table.read(obsid_src)
            q = Quaternion.Quat([ obs['ra_nom'], obs['dec_nom'], obs['roll_nom']])

            y, z = Ska.quatutil.radec2yagzag(evts['RA'], evts['DEC'], q)
            pos = dict(
                time=np.array(evts['time']),
                yag=np.array(y * 3600),
                zag=np.array(z * 3600))
            pos_pick = open(os.path.join(obsdir, 'released_pos.pkl'), 'w')
            cPickle.dump(pos, pos_pick)
            pos_pick.close()


    pos_files = glob("auto/obs*/released_pos.pkl")
    print "Retrieved sources for {} observations".format(len(pos_files))


GRADIENTS = dict(OOBAGRD3=dict(yag=6.98145650e-04,
                               zag=9.51578351e-05,
                               ),
                 OOBAGRD6=dict(yag=-1.67009240e-03,
                               zag=-2.79084775e-03,
                               ))



#for obsid in obsids:
#
#    print obsid
#    obs_info, src = find_obsid_src(obsid)
#
#
#    # position data
#    if not os.path.exists(os.path.join(obsdir, 'pos.pkl')) or REDO:
#        obs = obs_info
#        evts = Table.read(obsid_src)
#        import Quaternion
#        import Ska.quatutil
#        q = Quaternion.Quat([ obs['ra_nom'], obs['dec_nom'], obs['roll_nom']])
#
#        y, z = Ska.quatutil.radec2yagzag(evts['RA'], evts['DEC'], q)
#        pos = dict(
#            time=np.array(evts['time']),
#            yag=np.array(y * 3600),
#            zag=np.array(z * 3600))
#
#        pos_pick = open(os.path.join(obsdir, 'pos.pkl'), 'w')
#        cPickle.dump( pos, pos_pick)
#        pos_pick.close()
#
#        # retrieve gradient telemetry
#        tstart = evts['time'][0]
#        tstop = evts['time'][-1]
#        gradients = fetch.MSIDset(GRADIENTS.keys(), tstart-100, tstop+100)
#        for msid in gradients:
#            # filter bad telemetry in place
#            filter_bad_telem(gradients[msid])
#            times = gradients[msid].times
#            evt_idx = np.searchsorted(times, evts['time'])
#            # find a mean gradient, because this calibration is relative to mean 
#            mean_gradient = np.mean(gradients[msid].vals[evt_idx])
#            # and smooth the telemetry to deal with slow changes and large step sizes..
#            smooth_gradient = smooth(gradients[msid].vals)
#            y += (smooth_gradient[evt_idx] - mean_gradient) * GRADIENTS[msid]['yag']
#            z += (smooth_gradient[evt_idx] - mean_gradient) * GRADIENTS[msid]['zag']
#

