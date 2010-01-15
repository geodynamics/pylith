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

## @file playpen/powerlaw/powerlaw_gendb.py

## @brief Python application to compute power-law parameters used by
## PyLith, given a spatial database describing the temperature and
## another spatial database describing the laboratory-derived
## properties for the various materials. The output is another spatial
## database containing the power-law parameters for PyLith.

import numpy
import math

from pyre.applications.Script import Script as Application

class PowerLawApp(Application):
  """
  Python application to compute power-law parameters used by PyLith,
  given input spatial databases describing the temperature and the
  laboratory-derived properties for the various materials. The output is
  another spatial database containing the power-law parameters for PyLith.
  """
  
  ## \b Properties
  ## @li \b reference_value Indicates whether reference stress or 
  ##   reference strain rate is provided as input.
  ## @li \b reference_stress Value for reference stress.
  ## @li \b reference_strain_rate Value for reference strain rate.
  ##
  ## \b Facilities
  ## @li \b db_exponent Spatial db for power-law exponent, n.
  ## @li \b db_activation_energy Spatial db for activation energy, Q.
  ## @li \b db_temperature Spatial db for temperature, T.
  ## @li \b db_powerlaw_coefficient Spatial db for power-law coefficient, Ae.
  ## @li \b geometry Geometry for output database.
  ## @li \b iohandler Object for writing database.

  import pyre.inventory
  from pyre.units.pressure import MPa
  from pyre.units.time import s

  refSelection = pyre.inventory.str("reference_value",
                                    default="strain_rate",
                                    validator=pyre.inventory.choice(['stress',
                                                                     'strain_rate']))
  refSelection.meta['tip'] = "Indicates whether reference stress or " \
      "reference strain rate is provided as input."

  refStress = pyre.inventory.dimensional("reference_stress", default=1.0*MPa)
  refStress.meta['tip'] = "Reference stress value."
  
  refStrainRate = pyre.inventory.dimensional("reference_strain_rate",
                                             default=1.0e-6/s)
  refStrainRate.meta['tip'] = "Reference strain rate value."

  from spatialdata.spatialdb.SimpleDB import SimpleDB

  dbExponent = pyre.inventory.facility("db_exponent",
                                       family="spatial_database",
                                       factory=SimpleDB)
  dbExponent.meta['tip'] = "Spatial db for power-law exponent, n."

  dbActivationE = pyre.inventory.facility("db_activation_energy",
                                       family="spatial_database",
                                       factory=SimpleDB)
  dbActivationE.meta['tip'] = "Spatial db for activation energy, Q."

  dbTemperature = pyre.inventory.facility("db_temperature",
                                       family="spatial_database",
                                       factory=SimpleDB)
  dbTemperature.meta['tip'] = "Spatial db for temperature, T."

  dbAe = pyre.inventory.facility("db_powerlaw_coefficient",
                                 family="spatial_database",
                                 factory=SimpleDB)
  dbAe.meta['tip'] = "Spatial db for power-law coefficient, Ae."

  from spatialdata.spatialdb.generator.Geometry import Geometry
  geometry = pyre.inventory.facility("geometry", family="geometry",
                                     factory=Geometry)
  geometry.meta['tip'] = "Geometry for output database."

  from spatialdata.spatialdb.SimpleIOAscii import SimpleIOAscii
  iohandler = pyre.inventory.facility("iohandler", family="simpledb_io",
                                      factory=SimpleIOAscii)
  iohandler.meta['tip'] = "Object for writing database."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="powerlaw_gendb"):
    Application.__init__(self, name)

    return


  def main(self, *args, **kwds):
    """
    Application driver.
    """
    # Get output points
    self._info.log("Reading geometry.")
    self.geometry.read()
    points = self.geometry.vertices
    coordsys = self.geometry.coordsys

    (npoints, spaceDim) = points.shape
    refStrainRate = numpy.zeros( (npoints,), dtype=numpy.float64)
    refStress = numpy.zeros( (npoints,), dtype=numpy.float64)

    # Query databases to get inputs at output points
    self._info.log("Querying for parameters at output points.")
    n = self._queryDB(self.dbExponent, "power-law-exponent", points, coordsys)
    Q = self._queryDB(self.dbActivationE, "activation-energy", points, coordsys)
    #Ae = self._queryDB(self.dbAe, "power-law-coefficient", points, coordsys)
    Ae = self._queryDB(self.dbAe, "flow-constant", points, coordsys)
    T = self._queryDB(self.dbTemperature, "temperature", points, coordsys)

    # Compute power-law parameters
    self._info.log("Computing parameters at output points.")
    from pyre.handbook.constants.fundamental import R
    At = 3**(0.5*(n+1))/2.0 * Ae * numpy.exp(-Q/(R.value*T))
    
    if self.refSelection == "stress":
      refStress[:] = self.refStress.value
      refStrainRate = self.refStress.value**n / At
    elif self.refSelection == "strain_rate":
      refStrainRate[:] = self.refStrainRate.value
      refStress = (self.refStrainRate.value / At)**(1.0/n)
    else:
      raise ValueError("Invalid value (%s) for reference value." % \
                         self.refSelection)
    
    refStressInfo = {'name': "reference-stress",
                     'units': "Pa",
                     'data': refStress.flatten()}
    refStrainRateInfo = {'name': "reference-strain-rate",
                         'units': "1/s",
                         'data': refStrainRate.flatten()}
    exponentInfo = {'name': "powerlaw-exponent",
                    'units': "none",
                    'data': n.flatten()}

    # Write database
    self._info.log("Writing database.")
    data = {'points': points,
            'coordsys': coordsys,
            'data_dim': self.geometry.dataDim,
            'values': [refStressInfo, refStrainRateInfo, exponentInfo]}
    self.iohandler.write(data)
    return


  def _queryDB(self, db, valueName, points, cs):
    """
    Query spatial database 
    """

    (npoints, spaceDim) = points.shape
    data = numpy.zeros( (npoints,1), dtype=numpy.float64)
    err = numpy.zeros( (npoints,), dtype=numpy.int32)

    db.open()
    db.queryVals([valueName])
    db.multiquery(data, err, points, cs)
    db.close()
    errSum = numpy.sum(err)
    if errSum > 0:
      msg = "Query for %s failed at %d points.\n" \
          "Coordinates of points:\n" % (valueName, errSum)
      msg += "%s" % points[err,:]
      raise ValueError(msg)

    return data


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = PowerLawApp()
  app.run()

# End of file
