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

## @file unittests/libtests/feassemble/data/ElasticityExplicitApp.py

## @brief Python application for generating C++ data files for testing
## C++ elasticity integrator objects.

from ElasticityApp import ElasticityApp

import numpy
import feutils

# ----------------------------------------------------------------------

# ElasticityExplicitApp class
class ElasticityExplicitApp(ElasticityApp):
  """
  Python application for generating C++ data files for testing C++
  elasticity integrator objects.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="elasticityexplicitapp"):
    """
    Constructor.
    """
    ElasticityApp.__init__(self, name)

    self.normViscosity = 0.1
    return
  

  def main(self):
    """
    Run the application.
    """
    self._collectData()
    self._calculateResidual()
    self._calculateJacobian()
    self._calcDtStable()
    self._initData()
    self.data.write(self.name)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _calcDtStable(self):
    """
    Calculate stable time step for explicit time stepping.
    """
    vp = ((self.lameLambda + 2.0*self.lameMu) / self.density)**0.5
    minCellWidth = self.mesh.minCellWidth
    self.dtStableExplicit = minCellWidth / vp
    return


  def _calculateResidual(self):
    """
    Calculate contribution to residual of operator for integrator.
    """
    self.valsResidual = self.formulation.calculateResidual(self)
    if self.useGravity:
      self.valsResidual += self._calculateGravityVec()
    return


  def _calculateJacobian(self):
    """
    Calculate contribution to Jacobian matrix of operator for integrator.
    """
    self.valsJacobian = self.formulation.calculateJacobian(self)
    return


  def _initData(self):

    ElasticityApp._initData(self)
    
    # Calculated values
    self.data.addScalar(vtype="PylithScalar", name="_dtStableExplicit",
                       value=self.dtStableExplicit,
                       format="%16.8e");
    return


# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = ElasticityExplicitApp()
  app.run()


# End of file 
