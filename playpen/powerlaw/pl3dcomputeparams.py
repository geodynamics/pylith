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
## parameters for PyLith, as well as user-specified elastic properties. The
## generated spatial database will have the same dimensions and coordinate
## specifications as the input properties database.

import numpy

from pyre.applications.Script import Script as Application

class Pl3dComputeParams(Application):
  """
  Python application to compute power-law parameters used by PyLith,
  given input spatial databases describing the elastic properties, the
  temperature, and the laboratory-derived properties for the various
  materials. The output is another spatial database containing the power-law
  parameters for PyLith, as well as user-specified elastic properties. The
  generated spatial database will have the same dimensions and coordinate
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
    ## @li \b output_db_name File name for output spatial database.
    ##
    ## \b Facilities
    ## @li \b db_input_elastic_properties Database of input elastic properties.
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

    outputDbName = pyre.inventory.str("output_db_name",
                                      default="mat_powerlaw.spatialdb")
    outputDbName.meta['tip'] = "File name of output spatial database."
    
    dbTemplate = pyre.inventory.str(
      "db_template",
      default="db_input_power_law_properties",
      validator=pyre.inventory.choice(["db_input_power_law_properties",
                                       "db_input_elastic_properties",
                                       "db_input_temperature"]))
    dbTemplate.meta['tip'] = "Which input DB to use as a template for output."

    from spatialdata.spatialdb.SimpleDB import SimpleDB
    dbInputElasticProperties = pyre.inventory.facility(
      "db_input_elastic_properties",
      family="spatial_database",
      factory=SimpleDB)
    dbInputElasticProperties.meta['tip'] = "DB for input elastic properties."

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

    return


  def main(self):
    # import pdb
    # pdb.set_trace()
    self._getParameters()
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
    self.dbTemplate = self.inventory.dbTemplate
    self.outputDbName = self.inventory.outputDbName

    # Facilities
    self.dbInputElasticProperties = self.inventory.dbInputElasticProperties
    self.dbInputPowerLawProperties = self.inventory.dbInputPowerLawProperties
    self.dbInputTemperature = self.inventory.dbInputTemperature
    self.dbOutputProperties = self.inventory.dbOutputProperties

    return


  def _getParameters(self):
    """
    Get parameters from different databases and output power-law properties.
    """

    # Need to open all the input databases and then get the number of locations
    # from the template database.
    # We then loop over these locations and compute the power-law parameters
    # for each point.
    self.dbInputElasticProperties.open()
    self.dbInputPowerLawProperties.open()
    self.dbInputTemperature.open()
    temp_queryVals = ["temperature"]
    elastic_queryVals = ["density", "Vs", "Vp"]
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
