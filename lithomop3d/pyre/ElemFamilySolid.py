#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Charles A. Williams
#                       Rensselaer Polytechnic Institute
#             Copyright (C) 2006 Rensselaer Polytechnic Institute
#
# 
# 	Permission is hereby granted, free of charge, to any person
# 	obtaining a copy of this software and associated documentation
# 	files (the "Software"), to deal in the Software without
# 	restriction, including without limitation the rights to use,
# 	copy, modify, merge, publish, distribute, sublicense, and/or
# 	sell copies of the Software, and to permit persons to whom the
# 	Software is furnished to do so, subject to the following
# 	conditions:
# 
# 	The above copyright notice and this permission notice shall be
# 	included in all copies or substantial portions of the Software.
# 
# 	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# 	EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
# 	OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# 	NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
# 	HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# 	WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# 	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# 	OTHER DEALINGS IN THE SOFTWARE.
#         
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from ElemFamily import ElemFamily

import lithomop3d as lm3d

# ElemFamilySolid class
class ElemFamilySolid(ElemFamily):
  """Python manager for 3D solid elements."""

  # INVENTORY /////////////////////////////////////////////////////////////////

  class Inventory(ElemFamily.Inventory):
    """Python object for managing ElemFamilySolid facilities and properties."""

    import pyre.inventory
    

    label = pyre.inventory.str("label", default="")
    label.meta['tip'] = "Label for element family solid."

    # At present, it is not possible to mix different element types, but I
    # still think it makes sense to put integration info into the element
    # family.  Also, I plan to put back in the ability to mix elements.
    quadratureOrder = pyre.inventory.str("quadratureOrder", default="Full")
    quadratureOrder.validator = pyre.inventory.choice(["Full","Reduced","Selective"])
    quadratureOrder.meta['tip'] = "Integration order for elements in family."

    from MatIsoElastic import MatIsoElastic
    material = pyre.inventory.facility("material", factory=MatIsoElastic)
    material.meta['tip'] = "Material for a group of elements."
    
    # PUBLIC METHODS //////////////////////////////////////////////////////////

    def __init__(self, name="elemfamilysolid"):
      """Constructor."""
      Component.__init__(self, name, facility="elemfamilysolid")

      self.label = ""
      self.elemType = None
      self.quadratureOrder = ""
      self.material = {'materialType': "",
                       'numProps': 0,
                       'numState': 0,
                       'numState0': 0,
                       'fptrMatPrt': None,
                       'fptrElasMat': None,
                       'fptrElasStrs': None,
                       'fptrTdMatinit': None,
                       'fptrTdStrs': None,
                       'fptrTdStrsMat': None,
                       'fptrTdPrestrMat': None,
                       'fptrGetState': None,
                       'fptrUpdateState': None}
      # Things that should be in ElemFamily:
      # ElemInfoSolid:  dimension and integration info
      # ElemFamilyState:stresses, strains, etc.
      # Materials:      material types (including all necessary function pointers),
      #                 database info?, property array.
      # Label:          descriptor for family
      # elemType:       string or integer description.
      # quadratureOrder:global for family
      
      self.elemFamily = None

      import journal
      self._info = journal.info(name)
      return

    # PRIVATE METHODS /////////////////////////////////////////////////////////

    def _configure(self):
      """Setup members using inventory."""
      self.label = self.inventory.label
      self.integrator = self.inventory.integrator
      self.material = self.inventory.material
      return

# version
__id__ = "$Id$"

# End of file 
