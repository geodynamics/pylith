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

from pylith.utils.CppDataNew import CppData

# ElasticityApp class
class ElasticityApp(CppData):
  """
  Python application for generating C++ data files for testing C++
  elasticity objects.
  """
  
  # Metadata
  dimension = 2
  
  # Test input data

  useInertia = False
  useBodyForce = False

  filenameMesh = "data/tri3.mesh"
  
  discretizationOrder = ["basisOrder", "quadOrder", "isContinuous"]
  discretizations = [
      {"basisOrder": 1, "quadOrder": 1, "isBasisContinuous": True}, # displacement
  ]
  
  lengthScale = 0
  timeScale = 0
  pressureScale = 0
  densityScale = 0


  def _collect(self):
    self.numSolnFields = len(self.discretizations)

    self._addScalar("char*", "_filenameMesh", self.filenameMesh, "\"%s\"")
    self._addScalar("bool", "_useInertia", self.useInertia, "%s")
    self._addScalar("bool", "_useBodyForce", self.useBodyForce, "%s")

    self._addScalar("int", "_numSolnFields", len(self.discretizations), "%d")
    #self._addStructArray("topology::Field::Discretization", "discretizations", self.discretizations, self.discretizationOrder)

    self._addScalar("PylithScalar", "_lengthScale", self.lengthScale, "%16.8e")
    self._addScalar("PylithScalar", "_timeScale", self.timeScale, "%16.8e")
    self._addScalar("PylithScalar", "_densityScale", self.densityScale, "%16.8e")
    self._addScalar("PylithScalar", "_pressureScale", self.pressureScale, "%16.8e")

    return


# ======================================================================
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("--namespace", action="store", dest="namespace", required=True)
    parser.add_argument("--objname", action="store", dest="objname", required=True)
    parser.add_argument("--parent", action="store", dest="parent")
    parser.add_argument("--header", action="store", dest="header", default="header.hh")

    args = parser.parse_args()

    app = ElasticityApp(args.namespace, args.objname, args.parent, args.header)
    app.creator = __file__
    app.run()

  
# End of file 
