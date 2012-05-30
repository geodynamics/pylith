import numpy
import h5py

# Load in change in tractions from coseismic simulation
filename = "output/step01-fault.h5"

h5 = h5py.File(filename, 'r', driver="sec2")
vertices = h5['geometry/vertices'][:]
tractions_change = h5['vertex_fields/traction_change'][0,:,:]
h5.close()

# Parameters for tractions associated with background stress field
# Nominal density
density = 2900.0
gacc = 9.80665
coef_friction = 0.6

# Background normal tractions are compressive (less than zero because
# y is negative)
tractions_bg_normal = density*gacc*(vertices[:,1])

# Background shear tractions are reverse ("right-lateral" is negative)
# because the normal tractions are negative.
tractions_bg_shear = coef_friction*tractions_bg_normal

tractions_shear = tractions_bg_shear + tractions_change[:,0]
tractions_normal = tractions_bg_normal + tractions_change[:,1]

filename = "afterslip_tractions.spatialdb"

# Create coordinate system for output
from spatialdata.geocoords.CSCart import CSCart
cs = CSCart()
cs._configure()
cs.setSpaceDim(2)

# Create writer for spatialdata base file
from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
writer = SimpleIOAscii()
writer.inventory.filename = filename
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
