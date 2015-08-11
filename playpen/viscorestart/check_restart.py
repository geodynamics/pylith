#!/usr/bin/env python

material = "oceanmantle"

import numpy
import h5py

filename = "output/step05-%s.h5" % material
h5 = h5py.File(filename, "r")
stressA = h5['cell_fields/stress'][:]
strainA = h5['cell_fields/total_strain'][:]
strainViscousA = h5['cell_fields/viscous_strain'][:]
h5.close()

filename = "output/step06-%s.h5" % material
h5 = h5py.File(filename, "r")
stressB = h5['cell_fields/stress'][:]
strainB = h5['cell_fields/total_strain'][:]
strainViscousB = h5['cell_fields/viscous_strain'][:]
h5.close()

iB = 0
iA = iB+1
print stressA[iA,:,:]
print stressB[iB,:,:]
