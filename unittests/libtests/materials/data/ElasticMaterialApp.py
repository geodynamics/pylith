#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# <LicenseText>
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
    self.dimension = None
    self.numDBValues = None
    self.numParameters = None
    self.numParamValues = None
    self.dbValues = None
    self.parameterNames = None
    self.dbData = None
    self.parameterData = None

    # Elastic material information
    self.numLocs = None
    self.density = None
    self.strain = None
    self.stress = None
    self.elasticConsts = None
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
    self.data.addScalar(vtype="int", name="_dimension",
                        value=self.dimension,
                        format="%d")
    self.data.addScalar(vtype="int", name="_numDBValues",
                        value=self.numDBValues,
                        format="%d")
    self.data.addScalar(vtype="int", name="_numParameters",
                        value=self.numParameters,
                        format="%d")
    self.data.addArray(vtype="int", name="_numParamValues",
                        values=self.numParamValues,
                        format="%d", ncols=1)
    self.data.addArray(vtype="char*", name="_dbValues", values=self.dbValues,
                       format="\"%s\"", ncols=1)
    self.data.addArray(vtype="char*", name="_parameterNames",
                       values=self.parameterNames,
                       format="\"%s\"", ncols=1)
    self.data.addArray(vtype="double", name="_dbData", values=self.dbData,
                       format="%16.8e", ncols=1)
    self.data.addArray(vtype="double", name="_parameterData",
                       values=self.parameterData,
                       format="%16.8e", ncols=1)

    self.data.addScalar(vtype="int", name="_numLocs", value=self.numLocs,
                        format="%d")
    self.data.addArray(vtype="double", name="_density",
                       values=self.density,
                       format="%16.8e", ncols=1)
    self.data.addArray(vtype="double", name="_strain",
                       values=self.strain,
                       format="%16.8e", ncols=1)
    self.data.addArray(vtype="double", name="_stress",
                       values=self.stress,
                       format="%16.8e", ncols=1)
    self.data.addArray(vtype="double", name="_elasticConsts",
                       values=self.elasticConsts,
                       format="%16.8e", ncols=1)
      
    return

  
# End of file 
