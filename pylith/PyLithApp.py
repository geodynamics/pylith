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

## @file pylith/PyLithApp.py

## @brief Python PyLith application

#from mpi.Application import Application
from pyre.applications.Script import Script as Application

# PyLithApp class
class PyLithApp(Application):
  """
  Python PyLithApp application.
  """
  
  # INVENTORY //////////////////////////////////////////////////////////

  class Inventory(Application.Inventory):
    """
    Python object for managing PyLithApp facilities and properties.
    """

    ## @class Inventory
    ## Python object for managing PyLithApp facilities and properties.
    ##
    ## \b Properties
    ## @li None
    ##
    ## \b Facilities
    ## @li \b mesher Generates or imports the computational mesh.
    ## @li \b problem Computational problem to solve

    import pyre.inventory

    from pyre.units.time import second
    totalTime = pyre.inventory.dimensional("total_time", default=0.0*second,
                          validator=pyre.inventory.greaterEqual(0.0*second))
    totalTime.meta['tip'] = "Time duration for simulation."

    from pylith.topology.MeshImporter import MeshImporter
    mesher = pyre.inventory.facility("mesh_generator", factory=MeshImporter)
    mesher.meta['tip'] = "Generates or imports the computational mesh."

    from pylith.problems.EqDeformation import EqDeformation
    problem = pyre.inventory.facility("problem", factory=EqDeformation)
    problem.meta['tip'] = "Computational problem to solve."


  # PUBLIC METHODS /////////////////////////////////////////////////////

  def main(self):
    """
    Run the application.
    """

    #mesh = self.mesher.create()
    #self.problem.mesh = mesh.distribute()
    self.problem.run(self)
    return
  

  def __init__(self, name="pylithapp"):
    """
    Constructor.
    """
    Application.__init__(self, name)
    self.mesher = None
    self.problem = None
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _configure(self):
    """
    Setup members using inventory.
    """
    Application._configure(self)
    self.mesher = self.inventory.mesher
    self.problem = self.inventory.problem
    return


# End of file 
