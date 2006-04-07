#!/bin/sh

lithomop3dapp.py -typos=relaxed -lm3dscan.fileRoot=lintet/bm5/bm5 -lm3dscan.ucdOutput=binary -log_summary -pc_type ilu -ksp_monitor -ksp_view -ksp_rtol 1e-9

# end of file
