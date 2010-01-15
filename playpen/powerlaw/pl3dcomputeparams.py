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

## @file powerlaw/pl3dcomputeparams

## @brief Python application to compute power-law parameters used by PyLith,
## given a spatial database describing the temperature and another spatial
## database describing the laboratory-derived properties for the various
## materials. The output is another spatial database containing the power-law
## parameters for PyLith. The generated spatial database will have the same
## dimensions and coordinate specifications as the input properties database.

import numpy
import math

from pyre.applications.Script import Script as Application

class Pl3dComputeParams(Application):
  """
  Python application to compute power-law parameters used by PyLith,
  given input spatial databases describing the temperature and the
  laboratory-derived properties for the various materials. The output is
  another spatial database containing the power-law parameters for PyLith.
  The generated spatial database will have the same dimensions and coordinate
  specifications as the input properties database.
  """
  
  class Inventory(Application.Inventory):
    """
    Python object for managing Pl3dComputeParams facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing Pl3dComputeParams facilities and properties.
    ##
    ## \b Properties
    ## @li \b use_reference_stress Use reference stress to compute reference strain rate.
    ## @li \b reference_stress Value of reference stress.
    ## @li \b reference_strain_rate Value of reference strain rate.
    ## @li \b db_template Which database to use as template for output.
    ## @li \b point_locations_file File containing points to output.
    ## @li \b output_db_name File name for output spatial database.
    ##
    ## \b Facilities
    ## @li \b db_input_power_law_properties Database of input power law properties.
    ## @li \b db_input_temperature Database of input temperature values.
    ## @li \b db_output_properties Database of output material properties.

    import pyre.inventory
    from pyre.units.pressure import Pa
    from pyre.units.time import s

    useReferenceStress = pyre.inventory.bool("use_reference_stress",
                                             default=False)
    useReferenceStress.meta['tip'] = "Use reference stress to compute reference strain rate."

    referenceStress = pyre.inventory.dimensional("reference_stress",
                                                 default=1.0*MPa)
    referenceStress.meta['tip'] = "Reference stress value."

    referenceStrainRate = pyre.inventory.dimensional("reference_strain_rate",
                                                     default=1.0e-6/s)
    referenceStrainRate.meta['tip'] = "Reference strain rate value."

    pointLocationsFile = pyre.inventory.str("point_locations_file",
                                            default="points.txt")
    pointLocationsFile.meta['tip'] = "File containing points to output."

    outputDbName = pyre.inventory.str("output_db_name",
                                      default="mat_powerlaw.spatialdb")
    outputDbName.meta['tip'] = "File name of output spatial database."
    
    from spatialdata.spatialdb.SimpleDB import SimpleDB

    dbInputPowerLawProperties = pyre.inventory.facility(
      "db_input_power_law_properties",
      family="spatial_database",
      factory=SimpleDB)
    dbInputPowerLawProperties.meta['tip'] = "DB for input power law properties."

    dbInputTemperature = pyre.inventory.facility(
      "db_input_temperature",
      family="spatial_database",
      factory=SimpleDB)
    dbInputTemperature.meta['tip'] = "DB for input temperature values."
    
    dbOutputProperties = pyre.inventory.facility(
      "db_output_properties",
      family="spatial_database",
      factory=SimpleDB)
    dbOutputProperties.meta['tip'] = "DB for output material properties."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="pl3dcomputeparams"):
    Application.__init__(self, name)

    self.numPoints = 0
    self.points = None

    return


  def main(self):
    # import pdb
    # pdb.set_trace()
    self._readPoints()
    self._computeParams()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Application._configure(self)

    # Parameters
    self.useReferenceStress = self.inventory.useReferenceStress
    self.referenceStress = self.inventory.referenceStress.value
    self.referenceStrainRate = self.inventory.referenceStrainRate.value
    self.pointLocationsFile = self.inventory.pointLocationsFile
    self.outputDbName = self.inventory.outputDbName

    # This script is presently restricted to 3D.
    self.spaceDim = 3

    # Facilities
    self.dbInputPowerLawProperties = self.inventory.dbInputPowerLawProperties
    self.dbInputTemperature = self.inventory.dbInputTemperature
    self.dbOutputProperties = self.inventory.dbOutputProperties

    return


  def _readPoints(self):
    """
    Read point coordinates to use in generating output.
    """
    inFile = file(self.pointLocationsFile, 'r')
    listCoords = []

    for line in inFile:
      if not line.strip():
        continue
      if re.compile('^#').search(line) is not None:
        continue
      else:
        self.numPoints += 1
        data = line.split()
        for dim in range(self.spaceDim):
          listCoords.append(float(data[dim]))

    self.points = numpy.array(
      listCoords, dtype=numpy.float64).reshape(self.numPoints, self.spaceDim) 

    return
    

  def _computeParams(self):
    """
    Get parameters from different databases and output power-law properties
    at specified points.
    """

    from pyre.handbook.constants.fundamental import R

    # Open input databases.
    self.dbInputPowerLawProperties.open()
    self.dbInputTemperature.open()

    temp_queryVals = ["temperature"]
    pl_queryVals = ["flow-constant", "activation-energy",
                    "power-law-exponent", "flow-constant-multiplier"]

    # Actually, I can't find a way of determining the number of locations
    # a priori.  We may need a different type of loop that continues until
    # the data in the template database is exhausted.
    # Alternatively, we can define a set of points for the output database
    # and just loop over those. One method might be:
    # specify min and max values for each coordinate direction along with
    # either an interval or number of increments. If the min and max values are
    # equal for a direction, that direction is just held constant.

    for loc in range(nlocs):
      
    return


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = Pl3dComputeParams()
  app.run()

# End of file
