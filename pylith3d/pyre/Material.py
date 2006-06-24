#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  PyLith by Charles A. Williams, Brad Aagaard, and Matt Knepley
#
#  Copyright (c) 2004-2006 Rensselaer Polytechnic Institute
#
#  Permission is hereby granted, free of charge, to any person obtaining
#  a copy of this software and associated documentation files (the
#  "Software"), to deal in the Software without restriction, including
#  without limitation the rights to use, copy, modify, merge, publish,
#  distribute, sublicense, and/or sell copies of the Software, and to
#  permit persons to whom the Software is furnished to do so, subject to
#  the following conditions:
#
#  The above copyright notice and this permission notice shall be
#  included in all copies or substantial portions of the Software.
#
#  THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
#  EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
#  MERCHANTABILITY,    FITNESS    FOR    A   PARTICULAR    PURPOSE    AND
#  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
#  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#  OF CONTRACT, TORT OR OTHERWISE,  ARISING FROM, OUT OF OR IN CONNECTION
#  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from pyre.components.Component import Component

# Material class
class Material(Component):
  """Python manager for material properties."""

  # INVENTORY //////////////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """Python object for managing Material facilities and properties."""

    import pyre.inventory

    from SpatialDB import SpatialDB
    db = pyre.inventory.facility("db", factory=SpatialDB)
    db.meta['tip'] = "Spatial database of material properties."

  # PUBLIC METHODS /////////////////////////////////////////////////////////////

  def getmodel(self):
    return self.materialModel

  def getmodeldb(self):
    return self.db
  
  def __init__(self, name="material"):
    """Constructor."""
    Component.__init__(self, name, facility="material")
    self.materialModel = {'materialType': "",
                          'numProps': 0,
                          'propNames': [""],
                          'numStateVars': 0,
                          'stateVarNames': [""],
                          'numState0Vars': 0,
                          'state0VarNames': [""],
                          'fptrMatPrt': None,
                          'fptrElasMat': None,
                          'fptrElasStrs': None,
                          'fptrTdMatinit': None,
                          'fptrTdStrs': None,
                          'fptrTdStrsMat': None,
                          'fptrTdPrestrMat': None,
                          'fptrGetState': None,
                          'fptrUpdateState': None}
    return

  def _configure(self):
    """Set configuration using inventory."""
    self.db = self.inventory.db
    return
  

# version
__id__ = "$Id$"

# End of file 
