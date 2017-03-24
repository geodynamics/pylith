#!/usr/bin/env python

material = "genmaxps"

import numpy
import h5py

filename = "output/grav_static_%s-visco.h5" % material
h5 = h5py.File(filename, "r")
stressA = h5['cell_fields/stress'][:]
strainA = h5['cell_fields/total_strain'][:]
strainViscous1A = h5['cell_fields/viscous_strain_1'][:]
strainViscous2A = h5['cell_fields/viscous_strain_2'][:]
strainViscous3A = h5['cell_fields/viscous_strain_3'][:]
h5.close()

filename = "output/grav_restart_%s-visco.h5" % material
h5 = h5py.File(filename, "r")
stressB = h5['cell_fields/stress'][:]
strainB = h5['cell_fields/total_strain'][:]
strainViscous1B = h5['cell_fields/viscous_strain_1'][:]
strainViscous2B = h5['cell_fields/viscous_strain_2'][:]
strainViscous3B = h5['cell_fields/viscous_strain_3'][:]
h5.close()

fields = {'stress': (stressA, stressB),
          'viscous_strain_1': (strainViscous1A, strainViscous1B),
          'viscous_strain_2': (strainViscous2A, strainViscous2B),
          'viscous_strain_3': (strainViscous3A, strainViscous3B) }

iB = 0
iA = iB+1
for name,fields in fields.items():
    fieldA,fieldB = fields
    small = 1.0e-6
    mask = numpy.abs(fieldA[iA,:,:]) > small
    norm = mask*fieldA[iA,:,:] + ~mask*small
    diff = (fieldA[iA,:,:] - fieldB[iB,:,:]) / norm
    if numpy.sum(numpy.abs(diff)) > 1.0e-4:
        print name
        print fieldA[iA,:,:],fieldB[iB,:,:],norm,diff
    else:
        print name,"OK"
