#!/usr/bin/env python
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
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file unittests/libtests/materials/data/IsotropicLinearElasticityPlaneStrain.py

## @brief Python application for generating C++ data files for testing
## C++ IsotropicLinearElasticityPlaneStrain objects.

import numpy

from pylith.tests.CppTestData import CppTestData

# ElasticityApp class
class ElasticityApp(CppTestData):
  """
  Python application for generating C++ data files for testing C++
  elasticity objects.
  """
  
  # Metadata
  dimension = 2
  cincludes = ["pylith/fekernels/elasticity.h"]
  
  # Test input data

  filenameMesh = "data/tri3.mesh"
  materialId = 0
  materialLabel = "IsotropicLinearElascitity"

  useInertia = False
  useBodyForce = False
  residualKernels = [
    "NULL", # f0
    "pylith_fekernels_f1_IsotropicLinearElasticityPlaneStrain" # f1
  ]
  jacobianKernels = [
    "NULL", # g0
    "NULL", # g1
    "NULL", # g2
    "pylith_fekernels_g3_uu_IsotropicLinearElasticityPlaneStrain", # g3
  ]
  
  lengthScale = 1000.0
  timeScale = 2.0
  pressureScale = 2.25e+10
  velScale = lengthScale / timeScale
  densityScale = pressureScale / velScale**2

  discretizations = [ # basisOrder, quadOrder, isBasisContinuous
      "{1, 1, true}", # displacement
  ]
  
  filenameAuxFieldsDB = "data/matinitialize.spatialdb"
  auxFields = [2500.0/densityScale, 2.25e+10/pressureScale, 2.25e+10/pressureScale]


  def __init__(self):
    CppTestData.__init__(self, "pylith.materials", "IsotropicLinearElasticityPlaneStrainData_Tri3", "IsotropicLinearElasticityPlaneStrainData")
    import os.path
    self.creator = os.path.relpath(__file__)
    return


  def _collect(self):
    self.numSolnFields = len(self.discretizations)

    self._addScalar("char*", "_filenameMesh", self.filenameMesh, "\"%s\"")
    self._addScalar("char*", "_label", self.materialLabel, "\"%s\"")
    self._addScalar("int", "_id", self.materialId, "%d")
    self._addScalar("int", "_dimension", self.dimension, "%d")

    self._addScalar("bool", "_useInertia", self.useInertia, "%s")
    self._addScalar("bool", "_useBodyForce", self.useBodyForce, "%s")

    self._addArray("PetscPointFunc", "_residualKernels", self.residualKernels, "  %s")
    self._addArray("PetscPointJac", "_jacobianKernels", self.jacobianKernels, "  %s")

    self._addScalar("char*", "_filenameAuxFieldsDB", self.filenameAuxFieldsDB, "\"%s\"")


    self._addScalar("int", "_numSolnFields", len(self.discretizations), "%d")
    self._addArray("pylith::topology::Field::DiscretizeInfo", "_discretizations", self.discretizations, "  %s")

    self._addArray("PylithScalar", "_auxFieldsFromDB", self.auxFields, "%16.8e")

    self._addScalar("PylithReal", "_lengthScale", self.lengthScale, "%16.8e")
    self._addScalar("PylithReal", "_timeScale", self.timeScale, "%16.8e")
    self._addScalar("PylithReal", "_densityScale", self.densityScale, "%16.8e")
    self._addScalar("PylithReal", "_pressureScale", self.pressureScale, "%16.8e")

    return


# ======================================================================
def generate():
  app = ElasticityApp()
  app.run()

  
# End of file 
