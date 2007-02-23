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

## @file pylith/materials/ElasticIsotropic3D.py

## @brief Python object for 3-D isotropic linear elastic material.

from Material import Material

# ElasticIsotropic3D class
class ElasticIsotropic3D(Material):
  """Python object for 3-D isotropic linear elastic constitutive model."""

  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Material.Inventory):
    """Python object for managing ElasticIsotropic3D facilities and
    properties."""

    ## @class Inventory
    ## Python object for managing ElasticIsotropic3D facilities and
    ## properties.
    ##
    ## \b Properties
    ## @li \b use_db Use spatial database for properties instead of
    ##               uniform values supplied here.
    ## @li \b vp P-wave speed.
    ## @li \b vs S-wave speed.
    ## @li \b density Mass density.
    ##
    ## \b Facilities
    ## @li None

    import pyre.inventory

    useDB = pyre.inventory.bool("use_db", default=False)
    useDB.meta['tip'] = "Use spatial database for properties instead of " \
                        "uniform values supplied here."

    from pyre.units.length import km
    from pyre.units.time import s
    vs = pyre.inventory.dimensional("vs", default=3.0*km/s)
    vs.meta['tip'] = "S-wave speed."

    vp = pyre.inventory.dimensional("vs", default=3.0*(3**0.5)*km/s)
    vp.meta['tip'] = "P-wave speed."

    from pyre.units.mass import kg
    from pyre.units.length import m
    density = pyre.inventory.dimensional("density", default=3000.0*kg/m**3)
    density.meta['tip'] = "Mass density."
    
  
  # PUBLIC METHODS /////////////////////////////////////////////////////

  def initialize(self):
    """
    Initialize material. If not using predefined database, create one
    with given parameters.
    """
    if not self.useDB:
      self._info.log("Creating trivial database for uniform properties.")
      # :TODO: Need method to create trivial database
      #from spatialdata.spatialdb.SimpleDB import createdb
      #import numpy
      #coords = numpy.zeros( (1,3), dtype=numpy.float64)
      #data = numpy.array( [ [self.vp, self.vs, self.density] ] )
      #self.db = createdb(values=["Vp", "Vs", "density"],
      #                   units=["m/s", "m/s", "kg/m^3"],
      #                   data=data,
      #                   coords=coords,
      #                   spaceDim=0)
    Material.initialize(self)
    return


  def __init__(self, name="elasticisotropic3d"):
    """
    Constructor.
    """
    Material.__init__(self, name)
    # :TODO: Need to create module for materials
    # import pylith.materials.materials as bindings
    # self.cppHandle = bindings.ElasticIsotropic3D()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """Set members using inventory."""
    self.useDB = self.inventory.useDB
    self.density = self.inventory.density
    self.vs = self.inventory.vs
    self.vp = self.inventory.vp
    return


# End of file 
