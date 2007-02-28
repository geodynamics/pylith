#!/bin/bash
# ======================================================================
#
# Shell script to run PyLith using small split node test.
#
# ======================================================================

# Do not use pythonTimestep for now until all the bugs are worked out.
#pyts="-pl3dscan.pythonTimestep=1"

set -x
pylith3dapp.py -pl3dscan.asciiOutput=full -pl3dscan.ucdOutput=ascii $pyts

# end of file
