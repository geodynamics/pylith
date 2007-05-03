#!/bin/sh

pylith3dapp.py -typos=relaxed -fileRoot=lintet/splittest/splittest -ucdOutput=binary -petsc.log_summary -petsc.pc_type=ilu -petsc.ksp_monitor -petsc.ksp_view -petsc.ksp_rtol=1e-9

# end of file
