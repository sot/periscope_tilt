"""Using ds9_contour.sh, make a good-sized thumbnail of each
   cut-out source for later visual inspection.  Expects
   sources in auto/obs* with region files "center.reg" and
   fits file "point_source.fits" for each obsid.  Makes
   a "ds9_src.png" in each obsid directory.
"""

from Ska.Shell import bash
import os
from glob import glob
import re

auto = glob("auto/obs*")
for dir in auto:
    if (os.path.exists(os.path.join(dir, 'center.reg')) and
        os.path.exists(os.path.join(dir, 'point_source.fits'))):
        if not os.path.exists(os.path.join(dir, 'ds9_src.png')):
            #print '%s/ds9_src.png' % dir
            #    bash('display obs%s/ds9_src.png' % ob) 
            regf = open(os.path.join(dir, 'center.reg'),'r')
            reg = regf.read()
            cmatch = re.match('circle\((\d+\.\d+),\s(\d+\.\d+),\s\d\)\n', reg)
            if cmatch:
                cenx = cmatch.group(1)
                ceny = cmatch.group(2)
                bash('./tilt/ds9_contour.sh %s/point_source.fits  %s %s %s/ds9_src.png'
                     % (dir, cenx, ceny, dir))
