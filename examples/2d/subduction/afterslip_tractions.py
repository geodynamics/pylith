# This python script extracts the changes in fault tractions from
# step01 and combines them with background depth-dependent fault
# tractions to create initial fault tractions to drive frictional
# afterslip.
#
# PREREQUISITES: numpy, h5py, output from step01
#
# Run the script:
#  python afterslip_tractions.py

import numpy
import h5py

# Load in change in tractions from coseismic simulation
h5 = h5py.File("output/step01-fault.h5", 'r', driver="sec2")
vertices = h5['geometry/vertices'][:]
tractions_change = h5['vertex_fields/traction_change'][0,:,:]
h5.close()

# Parameters for tractions associated with background stress field
# Nominal density
density = 2900.0
gacc = 9.80665
coef_friction = 0.6

# Background normal tractions are overburden and compressive
# (negative, y is negative)
tractions_bg_normal = density*gacc*(vertices[:,1])

# Background shear tractions are reverse (in 2-D right-lateral is negative)
# because the normal tractions are negative.
tractions_bg_shear = coef_friction*tractions_bg_normal

# Combine traction changes and background tractions
tractions_shear = tractions_bg_shear + tractions_change[:,0]
tractions_normal = tractions_bg_normal + tractions_change[:,1]

# Create coordinate system for spatial database
from spatialdata.geocoords.CSCart import CSCart
cs = CSCart()
cs._configure()
cs.setSpaceDim(2)

# Create writer for spatial database file
from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
writer = SimpleIOAscii()
writer.inventory.filename = "afterslip_tractions.spatialdb"
writer._configure()
writer.write({'points': vertices,
              'coordsys': cs,
              'data_dim': 1,
              'values': [{'name': "traction-shear",
                          'units': "Pa",
                          'data': tractions_shear},
                         {'name': "traction-normal",
                          'units': "Pa",
                          'data': tractions_normal}]})

# End of file
