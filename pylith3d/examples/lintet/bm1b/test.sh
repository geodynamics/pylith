#!/bin/sh

pylith3dapp.py -typos=relaxed -pl3dscan.fileRoot=lintet/bm1b/bm1b -log_summary -petsc_solver 1 -pc_type ilu -ksp_monitor -ksp_view -ksp_rtol 1e-9

# end of file
