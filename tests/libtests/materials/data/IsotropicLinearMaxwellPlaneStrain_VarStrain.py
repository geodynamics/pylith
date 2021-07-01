#!/usr/bin/env nemesis
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University at Buffalo
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2021 University of California, Davis
#
# See LICENSE.md for license information.
#
# ----------------------------------------------------------------------
#

# @file tests/libtests/materials/data/IsotropicLinearMaxwellPlaneStrain_VarStrain.py

# @brief Python application for generating spatial database files for
# testing IsotropicLinearMaxwellPlaneStrain via Method of
# Manufactured Solutions for linearly varying total strain.

# ----------------------------------------------------------------------

# Domain
from spatialdata.geocoords.CSCart import CSCart
from spatialdata.spatialdb.SimpleGridAscii import createWriter
import math
import numpy
XLIM = (-4.0e+3, +4.0e+3)
YLIM = XLIM
DX = 100.0

# Material properties
DENSITY = 4000.0
VS = 5600.0
VP = 10000.0
VISCOSITY = 7.91700159488e+19

# Variable definitions for solution.
A = 1.0e-6
B = 2.5e-6
C = 3.0e-6
TIME = 1.0e7

# ----------------------------------------------------------------------


x = numpy.arange(XLIM[0], XLIM[1] + 0.1 * DX, DX, dtype=numpy.float64)
y = numpy.arange(YLIM[0], YLIM[1] + 0.1 * DX, DX, dtype=numpy.float64)
xgrid, ygrid = numpy.meshgrid(x, y)
points = numpy.vstack((xgrid.ravel(), ygrid.ravel())).transpose()
npts = points.shape[0]
PX = points[:, 0]
PY = points[:, 1]

density = DENSITY * numpy.ones((npts,))
vs = VS * numpy.ones((npts,))
vp = VP * numpy.ones((npts,))
viscosity = VISCOSITY * numpy.ones((npts,))

# Create material properties for solution.
shearModulus = DENSITY * VS * VS
lameConstant = DENSITY * VP * VP - 2.0 * shearModulus
bulkModulus = lameConstant + 2.0 * shearModulus / 3.0
maxwellTime = VISCOSITY / shearModulus

# Create coordinate system for spatial database
cs = CSCart()
cs._configure()
cs.setSpaceDim(2)

# ----------------------------------------------------------------------


def generateAuxSubfields():
    totalStrain_11 = (2.0 * A * PX + B * PY) * math.exp(-TIME / maxwellTime)
    totalStrain_12 = (B * PX / 2.0 + B * PY / 2.0 + C * PX +
                      C * PY) * math.exp(-TIME / maxwellTime)
    totalStrain_22 = (2.0 * A * PY + B * PX) * math.exp(-TIME / maxwellTime)
    totalStrain_33 = numpy.zeros_like(totalStrain_22)

    visStrain_11 = (math.exp(TIME / maxwellTime) - 1.0) * (A * PX - A * PY - B * PX + B * PY) * \
        math.exp(-2.0 * TIME / maxwellTime)
    visStrain_12 = (math.exp(TIME / maxwellTime) - 1.0) * (B * PX + B * PY + C * PX + C * PY) * \
        math.exp(-2.0 * TIME / maxwellTime)
    visStrain_22 = -(math.exp(TIME / maxwellTime) - 1.0) * (A * PX - A * PY - B * PX + B * PY) * \
        math.exp(-2.0 * TIME / maxwellTime)
    visStrain_33 = -(math.exp(TIME / maxwellTime) - 1.0) * (A * PX + A * PY + B * PX + B * PY) * \
        math.exp(-2.0 * TIME / maxwellTime)

    equil_1 = 4.0 * bulkModulus * \
        (A + B) * math.exp(-TIME / maxwellTime) * \
        numpy.ones(npts, dtype=numpy.float64)
    equil_2 = 4.0 * bulkModulus * \
        (A + B) * math.exp(-TIME / maxwellTime) * \
        numpy.ones(npts, dtype=numpy.float64)

    writer = createWriter(
        "IsotropicLinearMaxwellPlaneStrain_VarStrain_aux.spatialdb")
    writer.write({'points': points,
                  'x': x,
                  'y': y,
                  'coordsys': cs,
                  'data_dim': 2,
                  'values': [{'name': "vs", 'units': "m/s", 'data': vs},
                             {'name': "vp", 'units': "m/s", 'data': vp},
                             {'name': "density", 'units': "kg/m**3", 'data': density},
                             {'name': "viscosity", 'units': "Pa*s",
                                 'data': viscosity},
                             {'name': "total_strain_xx", 'units': "None",
                                 'data': totalStrain_11},
                             {'name': "total_strain_yy", 'units': "None",
                                 'data': totalStrain_22},
                             {'name': "total_strain_zz", 'units': "None",
                                 'data': totalStrain_33},
                             {'name': "total_strain_xy", 'units': "None",
                                 'data': totalStrain_12},
                             {'name': "vis_strain_xx", 'units': "None",
                                 'data': visStrain_11},
                             {'name': "vis_strain_yy", 'units': "None",
                                 'data': visStrain_22},
                             {'name': "vis_strain_zz", 'units': "None",
                                 'data': visStrain_33},
                             {'name': "vis_strain_xy", 'units': "None",
                                 'data': visStrain_12},
                             {'name': "body_force_x", 'units': "N", 'data': equil_1},
                             {'name': "body_force_y", 'units': "N", 'data': equil_2},
                             ]})
    return


# ----------------------------------------------------------------------
def generateSolution():

    disp = numpy.zeros((npts, 2))
    disp[:, 0] = (A * PX**2 + 2.0 * B * PX * PY + C * PY**2) * \
        math.exp(-TIME / maxwellTime)
    disp[:, 1] = (A * PY**2 + 2.0 * B * PX * PY + C * PX**2) * \
        math.exp(-TIME / maxwellTime)

    disp_dot = numpy.zeros((npts, 2))
    disp_dot[:, 0] = -(A * PX**2 + 2.0 * B * PX * PY + C *
                       PY**2) * math.exp(-TIME / maxwellTime) / maxwellTime
    disp_dot[:, 1] = -(A * PY**2 + 2.0 * B * PX * PY + C *
                       PX**2) * math.exp(-TIME / maxwellTime) / maxwellTime

    # Create writer for spatial database file
    writer = createWriter(
        "IsotropicLinearMaxwellPlaneStrain_VarStrain_soln.spatialdb")
    writer.write({'points': points,
                  'x': x,
                  'y': y,
                  'coordsys': cs,
                  'data_dim': 2,
                  'values': [{'name': "displacement_x", 'units': "m", 'data': disp[:, 0]},
                             {'name': "displacement_y",
                                 'units': "m", 'data': disp[:, 1]},
                             {'name': "displacement_dot_x",
                                 'units': "m/s", 'data': disp_dot[:, 0]},
                             {'name': "displacement_dot_y",
                                 'units': "m/s", 'data': disp_dot[:, 1]},
                             ]})

    return


# ----------------------------------------------------------------------
def generatePerturbation():
    PERT_DX = 500.0
    PERT_AMPLITUDE = 1.0e-2

    x = numpy.arange(XLIM[0], XLIM[1] + 0.1 * PERT_DX,
                     PERT_DX, dtype=numpy.float64)
    y = numpy.arange(YLIM[0], YLIM[1] + 0.1 * PERT_DX,
                     PERT_DX, dtype=numpy.float64)
    xgrid, ygrid = numpy.meshgrid(x, y)
    points = numpy.vstack((xgrid.ravel(), ygrid.ravel())).transpose()
    npts = points.shape[0]

    disp = PERT_AMPLITUDE * (numpy.random.rand(npts, 2) - 0.5)
    disp_dot = 0 * disp

    # Create writer for spatial database file
    writer = createWriter(
        "IsotropicLinearMaxwellPlaneStrain_VarStrain_pert.spatialdb")
    writer.write({'points': points,
                  'x': x,
                  'y': y,
                  'coordsys': cs,
                  'data_dim': 2,
                  'values': [{'name': "displacement_x", 'units': "m", 'data': disp[:, 0]},
                             {'name': "displacement_y",
                                 'units': "m", 'data': disp[:, 1]},
                             {'name': "displacement_dot_x",
                                 'units': "m/s", 'data': disp_dot[:, 0]},
                             {'name': "displacement_dot_y",
                                 'units': "m/s", 'data': disp_dot[:, 1]},
                             ]})

    return


# ======================================================================
def generate():
    generateAuxSubfields()
    generateSolution()
    generatePerturbation()


# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":
    generate()

# End of file
