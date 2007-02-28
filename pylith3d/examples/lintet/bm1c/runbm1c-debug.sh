#!/bin/bash
# ======================================================================
#
# Shell script to run PyLith using SCEC benchmark 1 (low resolution,
# 2 materials).
#
# ======================================================================

set -x
pylith3dapp.py -pl3dscan.asciiOutput=full -pl3dscan.ucdOutput=ascii $pyts

# end of file
