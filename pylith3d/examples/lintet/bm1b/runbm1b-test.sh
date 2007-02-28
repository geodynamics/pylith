#!/bin/bash
# ======================================================================
#
# Shell script to run PyLith using SCEC benchmark 1 (low resolution).
#
# ======================================================================

# Do not use pythonTimestep for now until all the bugs are worked out.
#pyts="-pl3dscan.pythonTimestep=1"

set -x
pylith3dapp.py -pl3dscan.asciiOutput=echo -pl3dscan.ucdOutput=ascii $pyts

# end of file
