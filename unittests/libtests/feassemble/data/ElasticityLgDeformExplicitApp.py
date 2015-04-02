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

## @file unittests/libtests/feassemble/data/ElasticityLgDeformExplicitApp.py

## @brief Python application for generating C++ data files for testing
## C++ elasticity integrator objects.

from ElasticityLgDeformApp import ElasticityLgDeformApp
from ElasticityExplicitApp import ElasticityExplicitApp

import numpy

# ----------------------------------------------------------------------

# ElasticityLgDeformExplicitApp class
class ElasticityLgDeformExplicitApp(ElasticityLgDeformApp, 
                                    ElasticityExplicitApp):
  """
  Python application for generating C++ data files for testing C++
  elasticity integrator objects.
  """
  
  # PUBLIC METHODS /////////////////////////////////////////////////////
  
  def __init__(self, name="elasticitylgdeformexplicitapp"):
    """
    Constructor.
    """
    ElasticityLgDeformApp.__init__(self, name)
    ElasticityExplicitApp.__init__(self, name)
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
  

# MAIN /////////////////////////////////////////////////////////////////
if __name__ == "__main__":

  app = ElasticityLgDeformExplicitApp()
  app.run()


# End of file 
