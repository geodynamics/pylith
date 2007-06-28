#!/bin/bash
pylith axialdisp.cfg 2>&1 | tee axialdisp-np1.log
mv axialdisp_t0.vtk axialdisp_t0.np1.vtk
pylith axialdisp.cfg --nodes=2 2>&1 | tee axialdisp-np2.log
mv axialdisp_t0.vtk axialdisp_t0.np2.vtk
pylith dislocation.cfg 2>&1 | tee dislocation-np1.log
mv dislocation_t0.vtk dislocation_t0.np1.vtk
pylith dislocation.cfg --nodes=2 2>&1 | tee dislocation-np2.log
mv dislocation_t0.vtk dislocation_t0.np2.vtk
pylith sheardisp.cfg 2>&1 | tee sheardisp-np1.log
mv sheardisp_t0.vtk sheardisp_t0.np1.vtk
pylith sheardisp.cfg --nodes=2 2>&1 | tee sheardisp-np2.log
mv sheardisp_t0.vtk sheardisp_t0.np2.vtk
