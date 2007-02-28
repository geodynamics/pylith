#!/bin/sh

pylith3dapp.py -pl3dscan.fileRoot=lintet/bm5/bm5 -pl3dscan.ucdOutput=binary -log_summary -pc_type ilu -ksp_monitor -ksp_view -ksp_rtol 1e-9

# end of file
