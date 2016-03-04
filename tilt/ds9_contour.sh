#!/bin/bash

ds9 -title ds9 \
-file $1 \
-bin about $2 $3 \
-bin factor .2 \
-zoom 4 \
-scale sqrt \
-frame new \
-file $1 \
-bin about $2 $3 \
-scale sqrt \
-tile yes \
-contour nlevels 10 \
-contour yes \
-saveimage png $4 \
-exit
