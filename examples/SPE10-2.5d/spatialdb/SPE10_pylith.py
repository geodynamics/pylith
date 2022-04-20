#!/usr/bin/env nemesis
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2016 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------

import numpy 
from netCDF4 import Dataset
import scipy.ndimage
import matplotlib.pyplot as plt
from spatialdata.geocoords.CSCart import CSCart

# mesh = Dataset('mesh_hex.exo','r')

# SPE 10 Parameters
nx = 60
ny = 220
nz = 85

dx = 20.
dy = 10. 
dz = 2. 

# Select which layer in xy plane
zval = 0

phi = numpy.loadtxt('spe_phi.dat').ravel()
phi = phi.reshape([nz,ny,nx])

perm = numpy.loadtxt('spe_perm.dat').ravel() * 9.869233e-16 

k_x = perm[0: nz*ny*nx].reshape([nz,ny,nx])
k_y = perm[nz*ny*nx:2*nz*ny*nx].reshape([nz,ny,nx])
k_z = perm[2*nz*ny*nx:].reshape([nz,ny,nx])

k_x = k_x[zval,:,:]
k_y = k_y[zval,:,:]
k_z = k_z[zval,:,:]

k_tensor = numpy.zeros([nx*ny,4])
k_tensor[:,0] = k_x.flatten()
k_tensor[:,1] = k_y.flatten()
k_tensor[:,2] = k_z.flatten()
k_tensor[:,3] = numpy.zeros(nx*ny)
k_tensor[:,4] = numpy.zeros(nx*ny)
k_tensor[:,5] = numpy.zeros(nx*ny)

phi = phi[zval,:,:]

# Two Dimensional Values

fluid_density = 1000 * numpy.ones(nx*ny) # kg / m**3
solid_density = 2300 * numpy.ones(nx*ny)# kg / m**3
biot_coefficient = 0.771 * numpy.ones(nx*ny)
fluid_viscosity = 0.001 * numpy.ones(nx*ny) # Pa*s
shear_modulus = 8662741799.83 * numpy.ones(nx*ny) # Pa
fluid_bulk_modulus = 2e9 * numpy.ones(nx*ny) # Pa
drained_bulk_modulus = 10e9 * numpy.ones(nx*ny) # Pa
solid_bulk_modulus = 11039657020.4 * numpy.ones(nx*ny) # Pa

class GenerateDB(object):
  def run(self):
    """Generate the database.
    """
        
    # Domain, centroid
    x = (numpy.arange(0, nx*dx, dx) + dx/2)
    y = (numpy.arange(0, ny*dy, dy) + dy/2)
    z = (numpy.arange(0, nz*dz, dz) + dz/2)
    
    x3, y3, z3 = numpy.meshgrid(x,y,z)
    
    npts_x = x.shape[0]
    npts_y = y.shape[0]
    npts_z = z.shape[0]
    
    xyz = numpy.zeros((npts_x*npts_y*npts_z, 3), dtype=numpy.float64)
    xyz[:, 0] = x3.ravel()
    xyz[:, 1] = y3.ravel()
    xyz[:, 2] = z3.ravel()
    
    cs = CSCart()
    cs.inventory.units = 'foot'
    cs.inventory.spaceDim = 3
    cs._configure()
    
    porosityGrid = {"name": "porosity",
                    "units": "none",
                    "data": numpy.ravel(phi)}

    permeability_xxGrid = {"name": "permeability_xx",
                           "units": "m*m",
                           "data": k_tensor[:,0]}        

    permeability_yyGrid = {"name": "permeability_yy",
                           "units": "m*m",
                           "data": k_tensor[:,1]}

    permeability_zzGrid = {"name": "permeability_zz",
                           "units": "m*m",
                           "data": k_tensor[:,2]}

    permeability_xyGrid = {"name": "permeability_xy",
                           "units": "m*m",
                           "data": k_tensor[:,3]}
                           
    permeability_yzGrid = {"name": "permeability_yz",
                           "units": "m*m",
                           "data": k_tensor[:,4]}
                           
    permeability_xzGrid = {"name": "permeability_xz",
                           "units": "m*m",
                           "data": k_tensor[:,5]}                               
                           
    solid_densityGrid = {"name": "solid_density",
                         "units": "kg/m**3",
                         "data": solid_density.flatten()}                    

    fluid_densityGrid = {"name": "fluid_density",
                         "units": "kg/m**3",
                         "data": fluid_density.flatten()}

    fluid_viscosityGrid = {"name": "fluid_viscosity",
                           "units": "Pa*s",
                           "data": fluid_viscosity.flatten()}

    shear_modulusGrid = {"name": "shear_modulus",
                         "units": "Pa",
                         "data": shear_modulus.flatten()}
  
    drained_bulk_modulusGrid = {"name": "drained_bulk_modulus",
                                "units": "Pa",
                                "data": drained_bulk_modulus.flatten()}
                                
    solid_bulk_modulusGrid = {"name": "solid_bulk_modulus",
                              "units": "Pa",
                              "data": solid_bulk_modulus.flatten()}

    fluid_bulk_modulusGrid = {"name": "fluid_bulk_modulus",
                              "units": "Pa",                 
                              "data": fluid_bulk_modulus.flatten()}
    
    biot_coefficientGrid = {"name": "biot_coefficient",
                            "units": "none",
                            "data": biot_coefficient.flatten()}

    data = {"num-x": x.shape[0],
            "num-y": y.shape[0]
            "num-z": z.shape[0],
            "points": xyz,                
            "x": x,
            "y": y,
            "z": z,
            "coordsys": cs,
            "data_dim": 3,
            "values": [porosityGrid, permeability_xxGrid, permeability_yyGrid,
            permeability_zzGrid, permeability_xyGrid, permeability_yzGrid, 
            permeability_xzGrid, solid_densityGrid, fluid_densityGrid,
            fluid_viscosityGrid, shear_modulusGrid, drained_bulk_modulusGrid,
            solid_bulk_modulusGrid, fluid_bulk_modulusGrid, biot_coefficientGrid]}

    from spatialdata.spatialdb.SimpleIOAscii import createWriter
    io = createWriter("SPE10_parameters.spatialdb")
    io.write(data)        
    return

# ======================================================================
if __name__ == "__main__":
    GenerateDB().run()


# End of file
