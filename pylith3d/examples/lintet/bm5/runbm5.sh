#!/bin/bash
# ======================================================================
#
# Shell script to run PyLith using SCEC benchmark 5 (low resolution).
#
# ======================================================================

# Do not use pythonTimestep for now until all the bugs are worked out.
#pyts="-pythonTimestep=1"

set -x
pylith3dapp.py -asciiOutput=none -ucdOutput=binary $pyts

# end of file
