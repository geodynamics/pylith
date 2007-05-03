#!/bin/bash
# ======================================================================
#
# Shell script to run PyLith using SCEC benchmark 1 (low resolution,
# 2 materials).
#
# ======================================================================

# Do not use pythonTimestep for now until all the bugs are worked out.
#pyts="-pythonTimestep=1"

set -x
pylith3dapp.py -asciiOutput=none -ucdOutput=binary $pyts

# end of file
