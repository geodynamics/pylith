#!/bin/bash
# ======================================================================
#
# Shell script to run PyLith using SCEC benchmark 1 (low resolution,
# 2 materials).
#
# ======================================================================

set -x
pylith3dapp.py -asciiOutput=full -ucdOutput=ascii $pyts

# end of file
