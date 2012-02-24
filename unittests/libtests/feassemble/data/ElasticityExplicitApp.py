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
# Copyright (c) 2010-2012 University of California, Davis
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
    self._calculateResidualLumped()
    self._calculateJacobianLumped()
    self._initData()
    self.data.write(self.name)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _calculateResidualLumped(self):
    """
    Calculate contribution to residual of operator for integrator.
    """
    self.valsResidualLumped = self.formulation.calculateResidualLumped(self)
    if self.useGravity:
      self.valsResidualLumped += self._calculateGravityVec()
    return


  def _calculateJacobianLumped(self):
    """
    Calculate contribution to Jacobian matrix of operator for integrator.
    """
    self.valsJacobianLumped = self.formulation.calculateJacobianLumped(self)
    return


  def _initData(self):

    ElasticityApp._initData(self)
    # Calculated values
    self.data.addArray(vtype="PylithScalar", name="_valsResidualLumped",
                       values=self.valsResidualLumped,
                       format="%16.8e", ncols=self.spaceDim)
    self.data.addArray(vtype="PylithScalar", name="_valsJacobianLumped",
                       values=self.valsJacobianLumped,
                       format="%16.8e", ncols=self.spaceDim)
    return


# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = ElasticityExplicitApp()
  app.run()


# End of file 
