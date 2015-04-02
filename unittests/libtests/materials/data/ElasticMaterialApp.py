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

## @file unittests/libtests/materials/data/ElasticMaterialApp.py

## @brief Python application for generating C++ data files for testing
## C++ elastic material objects.

from pyre.applications.Script import Script

import numpy

# ElasticMaterialApp class
class ElasticMaterialApp(Script):
  """
  Python application for generating C++ data files for testing C++
  elastic material objects.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Script.Inventory):
    """
    Python object for managing ElasticMaterialApp facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing ElasticMaterialApp facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b data Data manager.

    import pyre.inventory

    from pylith.utils.CppData import CppData
    data = pyre.inventory.facility("data", factory=CppData)
    data.meta['tip'] = "Data manager."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="elasticmaterialapp"):
    """
    Constructor.
    """
    Script.__init__(self, name)

    # Material information
    self.dimension = 0
    self.numLocs = 0
    self.numProperties = 0
    self.numStateVars = 0
    self.numDBProperties = 0
    self.numDBStateVars = 0
    self.numPropsQuadPt = 0
    self.numVarsQuadPt = 0
    self.numPropertyValues = None
    self.numStateVarValues = None
    self.dbPropertyValues = None
    self.dbStateVarValues = None
    self.dbProperties = None
    self.dbStateVars = None
    self.properties = None
    self.stateVars = None
    self.propertiesNondim = None
    self.stateVarsNondim = None
    self.lengthScale = 0
    self.timeScale = 0
    self.pressureScale = 0
    self.densityScale = 0

    # Elastic material information
    self.dtStableImplicit = 1.0e+99
    self.dtStableExplicit = 1.0e+99
    self.density = None
    self.strain = None
    self.stress = None
    self.elasticConsts = None
    self.initialStress = None
    self.initialStrain = None
    self.stateVarsUpdated = None
    return


  def main(self):
    """
    Run the application.
    """
    self._initData()
    self.data.write(self.name)
    return
  

  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Set members using inventory.
    """
    Script._configure(self)
    self.data = self.inventory.data
    return


  def _initData(self):
    self.numDBProperties = len(self.dbPropertyValues)
    if not self.dbStateVarValues is None:
      self.numDBStateVars = len(self.dbStateVarValues)
    self.numPropsQuadPt = numpy.sum(self.numPropertyValues)
    if not self.numStateVarValues is None:
      self.numVarsQuadPt = numpy.sum(self.numStateVarValues)
    self.numProperties = self.numPropertyValues.shape[0]
    if not self.numStateVarValues is None:
      self.numStateVars = self.numStateVarValues.shape[0]

    self.data.addScalar(vtype="int", name="_dimension",
                        value=self.dimension,
                        format="%d")
    self.data.addScalar(vtype="int", name="_numLocs",
                        value=self.numLocs,
                        format="%d")
    self.data.addScalar(vtype="int", name="_numProperties",
                        value=self.numProperties,
                        format="%d")
    self.data.addScalar(vtype="int", name="_numStateVars",
                        value=self.numStateVars,
                        format="%d")
    self.data.addScalar(vtype="int", name="_numDBProperties",
                        value=self.numDBProperties,
                        format="%d")
    self.data.addScalar(vtype="int", name="_numDBStateVars",
                        value=self.numDBStateVars,
                        format="%d")
    self.data.addScalar(vtype="int", name="_numPropsQuadPt",
                        value=self.numPropsQuadPt,
                        format="%d")
    self.data.addScalar(vtype="int", name="_numVarsQuadPt",
                        value=self.numVarsQuadPt,
                        format="%d")
    self.data.addArray(vtype="int", name="_numPropertyValues",
                        values=self.numPropertyValues,
                        format="%d", ncols=1)
    self.data.addArray(vtype="int", name="_numStateVarValues",
                        values=self.numStateVarValues,
                        format="%d", ncols=1)
    self.data.addArray(vtype="char*", name="_dbPropertyValues",
                       values=self.dbPropertyValues,
                       format="\"%s\"", ncols=1)
    self.data.addArray(vtype="char*", name="_dbStateVarValues",
                       values=self.dbStateVarValues,
		       format="\"%s\"", ncols=1)
    self.data.addArray(vtype="PylithScalar", name="_dbProperties",
                       values=self.dbProperties,
                       format="%16.8e", ncols=1)
    self.data.addArray(vtype="PylithScalar", name="_dbStateVars",
                       values=self.dbStateVars,
		       format="%16.8e", ncols=1)
    self.data.addArray(vtype="PylithScalar", name="_properties",
                       values=self.properties,
                       format="%16.8e", ncols=1)
    self.data.addArray(vtype="PylithScalar", name="_stateVars",
                       values=self.stateVars,
                       format="%16.8e", ncols=1)
    self.data.addArray(vtype="PylithScalar", name="_propertiesNondim",
                       values=self.propertiesNondim,
                       format="%16.8e", ncols=1)
    self.data.addArray(vtype="PylithScalar", name="_stateVarsNondim",
                       values=self.stateVarsNondim,
                       format="%16.8e", ncols=1)
    self.data.addScalar(vtype="PylithScalar", name="_lengthScale",
                        value=self.lengthScale,
                        format="%16.8e")
    self.data.addScalar(vtype="PylithScalar", name="_timeScale",
                        value=self.timeScale,
                        format="%16.8e")
    self.data.addScalar(vtype="PylithScalar", name="_pressureScale",
                        value=self.pressureScale,
                        format="%16.8e")
    self.data.addScalar(vtype="PylithScalar", name="_densityScale",
                        value=self.densityScale,
                        format="%16.8e")

    self.data.addScalar(vtype="PylithScalar", name="_dtStableImplicit",
                        value=self.dtStableImplicit,
                        format="%16.8e")
    self.data.addScalar(vtype="PylithScalar", name="_dtStableExplicit",
                        value=self.dtStableExplicit,
                        format="%16.8e")
    self.data.addArray(vtype="PylithScalar", name="_density",
                       values=self.density,
                       format="%16.8e", ncols=1)
    self.data.addArray(vtype="PylithScalar", name="_strain",
                       values=self.strain,
                       format="%16.8e", ncols=1)
    self.data.addArray(vtype="PylithScalar", name="_stress",
                       values=self.stress,
                       format="%16.8e", ncols=1)
    self.data.addArray(vtype="PylithScalar", name="_elasticConsts",
                       values=self.elasticConsts,
                       format="%16.8e", ncols=1)
    self.data.addArray(vtype="PylithScalar", name="_initialStress",
                       values=self.initialStress,
                       format="%16.8e", ncols=1)
    self.data.addArray(vtype="PylithScalar", name="_initialStrain",
                       values=self.initialStrain,
                       format="%16.8e", ncols=1)
    self.data.addArray(vtype="PylithScalar", name="_stateVarsUpdated",
                       values=self.stateVarsUpdated,
                       format="%16.8e", ncols=1)
      
    return

  
# End of file 
