import os
from glob import glob
from Ska.DBI import DBI
from astropy.table import Table
from Ska.Shell import bash
import re
import termios, sys, os
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
        print 'display %s/ds9_src.png  ' % src_dir
        disp = bash('display -geometry 800x600+600+1000 %s/ds9_src.png &' % src_dir) 
        proc = re.match('\[\d+\]\s(\d+)', disp[0])
        print proc.group(1)
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
    for key in ['kalman_datestart', 'kalman_datestop', 'kalman_tstart', 'kalman_tstop', 'grating', 'readmode', 'datamode']:
        srcdict[key] = obs[0][key]
    srcdict['point_source'] = stat
    srcs.append(srcdict)

srcs = Table(srcs)
srcs.write("src_table.dat", format='ascii.tab')

