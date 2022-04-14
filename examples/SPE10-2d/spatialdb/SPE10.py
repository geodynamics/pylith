

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy
import scipy.ndimage

mesh = Dataset('../mesh/mesh_quad.exo','r')

# SPE 10 Parameters

nx = 60
ny = 220
nz = 85

dx = 20.
dy = 10. 
dz = 2. 

# Select which layer in xy plane
zval = 0


x_c = np.arange(0, nx*dx, dx) + dx/2
y_c = np.arange(0, ny*dy, dy) + dy/2
z_c = np.arange(0, nz*dz, dz) + dz/2

x_c_3, y_c_3, z_c_3 = np.meshgrid(x_c, y_c, z_c)

xyz_c = np.zeros((nx*ny*nz, 3), dtype=np.float64)
xyz_c[:, 0] = x_c_3.ravel()
xyz_c[:, 1] = y_c_3.ravel()
xyz_c[:, 2] = z_c_3.ravel()

x_n = np.arange(0, nx*dx + dx, dx)
y_n = np.arange(0, ny*dy + dy, dy)
z_n = np.arange(0, nz*dz + dy, dz)

x_n_3, y_n_3, z_n_3 = np.meshgrid(x_n, y_n, z_n)

phi = np.loadtxt('spe_phi.dat').ravel()
phi = phi.reshape([nz,ny,nx])

perm = np.loadtxt('spe_perm.dat').ravel() * 9.869233e-16 

k_x = perm[0: nz*ny*nx].reshape([nz,ny,nx])
k_y = perm[nz*ny*nx:2*nz*ny*nx].reshape([nz,ny,nx])

k_x = k_x[zval,:,:]
k_y = k_y[zval,:,:]
phi = phi[zval,:,:]


x_n = numpy.arange(0, nx*dx + dx, dx)
y_n = numpy.arange(0, ny*dy + dy, dy)

# Domain, centroid
x = numpy.arange(0, nx*dx, dx) + dx/2
y = numpy.arange(0, ny*dy, dy) + dy/2
x2, y2 = numpy.meshgrid(x,y)
npts_x = x.shape[0]
npts_y = y.shape[0]

# xx = x.reshape([1, x.shape[0]]) * numpy.ones([ny ,1])
# yy = y.reshape([y.shape[0], 1]) * numpy.ones([1, nx])
xy = numpy.zeros((npts_x*npts_y, 2), dtype=numpy.float64)
# xy[:, 0] = numpy.ravel(xx)
# xy[:, 1] = numpy.ravel(numpy.transpose(yy))
xy[:, 0] = x2.ravel()
xy[:, 1] = y2.ravel()


# Checkerboard

factor = 10
nx_check = np.int32(nx/factor)
ny_check = np.int32(ny/factor)

check = np.indices([ny_check,nx_check]).sum(axis=0) % 2
check = scipy.ndimage.zoom(check, factor, order=0)

