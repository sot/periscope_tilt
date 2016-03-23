import os
from glob import glob
from Ska.DBI import DBI
from astropy.table import Table
from Ska.Shell import bash
import re
import termios, sys, os
import numpy as np
import cPickle
import matplotlib.pyplot as plt

TERMIOS = termios
def getkey():
    fd = sys.stdin.fileno()
    old = termios.tcgetattr(fd)
    new = termios.tcgetattr(fd)
    new[3] = new[3] & ~TERMIOS.ICANON & ~TERMIOS.ECHO
    new[6][TERMIOS.VMIN] = 1
    new[6][TERMIOS.VTIME] = 0
    termios.tcsetattr(fd, TERMIOS.TCSANOW, new)
    c = None
    try:
        c = os.read(fd, 1)
    finally:
        termios.tcsetattr(fd, TERMIOS.TCSAFLUSH, old)
    return c

def yz_diff_sum(pos):
    y_hist = np.histogram(pos['yag'] - np.mean(pos['yag']),
                          bins=np.arange(-4, 4, .1),
                          density=True)
    z_hist = np.histogram(pos['zag'] - np.mean(pos['zag']),
                          bins=np.arange(-4, 4, .1),
                          density=True)
    return np.sum(np.abs(y_hist[0] - z_hist[0]))


def hist_center(yag, zag):
    h, y, z = np.histogram2d(yag,
                             zag,
                             bins=100)
    idx = divmod(h.argmax(), h.shape[1])
    return y[idx[0]], z[idx[1]]


def cart_to_polar(y, z):
    theta = np.sqrt(y**2 + z**2)
    phi = np.arctan2(z, y)
    return(theta, phi)

def radial_sym(y, z):
    yc, zc = hist_center(pos['yag'], pos['zag'])
    t, p = cart_to_polar(pos['yag'] - yc, pos['zag'] - zc)
    h = np.histogram(np.degrees(p), bins=np.arange(-180, 190, 30))
#    raise ValueError
    radtest =  (h[0] * 1.0 / len(p)) > (3 * 1.0 / len(h[0]))
    print h[0]
    centertest = (np.count_nonzero(t < 1.0) * 1.0 / len(t)) < .5
    center_hist = np.histogram(t, bins=np.arange(0, 5, .2))
    return np.any(radtest) | centertest | (np.argmax(center_hist[0]) == 0)

COUNT_LIMIT = 2500


acadb = DBI(server='sybase', dbi='sybase', user='aca_read')
DATADIR = 'auto'

obs_srcs = glob(os.path.join(DATADIR, "obs*/picked_src.dat"))
srcs = []
for src_file in obs_srcs:
    src = Table.read(src_file, format='ascii')
    if src['NET_COUNTS'] < COUNT_LIMIT:
        continue
    src_dir = os.path.dirname(src_file)
    stat_file = os.path.join(src_dir, 'point_stat.dat')
    stat = False
    if os.path.exists(stat_file):
        stat_text = open(stat_file).read().strip()
        if stat_text == 'True':
            stat = True
    else:
        if not os.path.exists("%s/ds9_src.png" % src_dir):
            continue
        print 'display %s/ds9_src.png  ' % src_dir
        disp = bash('display -geometry 800x600+600+1000 %s/ds9_src.png &' % src_dir) 
        proc = re.match('\[\d+\]\s(\d+)', disp[0])
        print proc.group(1)
        pos = cPickle.load(open("{}/pos.pkl".format(src_dir)))
        yc, zc = hist_center(pos['yag'] - np.mean(pos['yag']),
                             pos['zag'] - np.mean(pos['zag']))
        print "OFFSET is {}".format(yc**2 + zc**2)
        tf = radial_sym(pos['yag'], pos['zag'])
        print tf
        

        diff_sum = yz_diff_sum(pos)
        print "SUM is {}".format(diff_sum)


        c = getkey()
        if c == 'q':     ## break on a Return/Enter 
            break
        if c == 'n':
            bash('echo "False" > {}'.format(stat_file))
            print "updating %s NOT POINT" % src_dir
        if c == 'm':
            bash('echo "None" > {}'.format(stat_file))
            print "updating %s ..Eh " % src_dir
        if c == 'y':
            bash('echo "True" > {}'.format(stat_file))
            print "updating %s POINT" % src_dir
        os.kill(int(proc.group(1)), 1)


    src = Table.read(src_file, format='ascii.tab')
    # append via dictionaries to avoid dealing with fact that ascdsver is
    # read as a float on occasion
    srcdict = dict(zip(src[0].colnames, src[0].data))
    obs = acadb.fetchall("select * from observations where obsid = {}".format(srcdict['obsid']))
    for key in ['kalman_datestart', 'kalman_datestop', 'kalman_tstart', 'kalman_tstop', 'grating', 'readmode', 'datamode', 'detector', 'sim_z', 'sim_z_offset']:
        srcdict[key] = obs[0][key]
    tilt_db = acadb.fetchall(
        "select max_oobagrd3 - min_oobagrd3 as rd3diff, max_oobagrd6 - min_oobagrd6 as rd6diff from obs_periscope_tilt where obsid = {}".format(srcdict['obsid']))
    if not len(tilt_db):
        continue
    srcdict['oobagrd3_diff'] = tilt_db[0]['rd3diff']
    srcdict['oobagrd6_diff'] = tilt_db[0]['rd6diff']

    from Ska.Sun import pitch
    obs_pitch = pitch(obs['ra_pnt'], obs['dec_pnt'], obs['kalman_tstart'])
    srcdict['pitch'] = obs_pitch
    srcdict['point_source'] = stat
    srcs.append(srcdict)

srcs = Table(srcs)
srcs.write("src_table.dat", format='ascii.tab')

