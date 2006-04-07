#!/bin/sh

lithomop3dapp.py -typos=relaxed -lm3dscan.fileRoot=linhex/slip1/slip1 -lm3dscan.asciiOutput=full -lm3dscan.autoRotateSlipperyNodes=False -log_summary -petsc_solver 1 -pc_type ilu -ksp_monitor -ksp_view -ksp_rtol 1e-9

# end of file
