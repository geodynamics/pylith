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
## @brief Python object for 3-D isotropic linear elastic constitutive model.

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
    """Initialize material. If not using predefined database, create
    one with given parameters."""
    if not self.useDB:
      from spatialdata.spatialdb.SimpleDB import SimpleDB
      self.db = SimpleDB()
      self.db.create(values=["Vp", "Vs", "density"],
                     units=["m/s", "m/s", "kg/m^3"],
                     data=[ [self.vp, self.vs, self.density] ],
                     topology="point")
    Material.initialize(self)
    return


  def __init__(self, name="elasticisotropic3d"):
    """Constructor."""
    Material.__init__(self, name)
    # NEED TO CREATE HANDLE TO C++ BINDINGS
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """Set members using inventory."""
    self.useDB = self.inventory.useDB
    self.muLame = self.inventory.density*self.inventory.vs*self.inventory.vs
    self.lambdaLame = self.inventory.density*self.inventory.vp*self.inventory.vp - 2.0*self.muLame
    self.density = self.inventory.density
    return


 # version
__id__ = "$Id$"

# End of file 
