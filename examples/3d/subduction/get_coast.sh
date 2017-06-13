#!/bin/sh
# Script to get coastline from GMT and make a VTK file for visualization.
# Note that GMT4 is used. Parameters and Python script would need modification
# for GMT5.
# PROJ.4 is also required for coordinate conversions.

# Filenames.
txtwgs84="cascadia_wgs84.txt"
txttm="cascadia_tm.txt"
vtktm="cascadia_tm.vtk"

# Region.
region="-R-129/-117/41/50"

# Resolution: (f)ull, (h)igh, (i)ntermediate, (l)ow, or (c)rude.
reso="-Dh"

# Draw rivers option.
# rivers="-Ia"
rivers="-I1"

# Draw political boundaries option.
political="-Na"

# Get coastlines from pscoast.
pscoast $region $reso $rivers $political -Jx1.0d -W -m'#' > $txtwgs84

# Do coordinate conversion to local Transverse Mercator.
cs2cs +proj=lonlat +ellps=WGS84 +datum=WGS84 +to +proj=tmerc +lon_0=-122.6765 +lat_0=45.5231 +datum=WGS84 +k=0.9996 $txtwgs84 > $txttm

# Run seg2vtk.py script to convert pscoast info to VTK file.
./seg2vtk.py --in_file=$txttm --out_file=$vtktm --elevation=0.1

/bin/rm $txtwgs84 $txttm
