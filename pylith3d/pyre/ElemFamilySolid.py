#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#  PyLith by Charles A. Williams
#  Copyright (c) 2003-2006 Rensselaer Polytechnic Institute
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

from ElemFamily import ElemFamily

import pylith3d as pl3d

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

    def initialize(self, elements, family):
      """Set up storage and define family."""

      self._info.log("Initializing element family '%s'." \
                     self.label)

      # Get element family dimension and integration info.
      from ElemIntegrationSolid import ElemIntegrationSolid
      self.elemIntegrationSolid = ElemIntegrationSolid.getintegration(family["elemFamilyType"],
                                               quadratureOrder)
      
      # Transfer element node info from global elements into element family array using
      # FORTRAN routine.  The elements object is no longer needed after this has been done
      # for all families.
      pl3d.getfamily(elements['numElems'],
                     elements['elemNodeArraySize'],
                     elements['ptrElemNodeArray'],
                     elements['ptrElemTypes'],
                     elements['ptrElemIds'],
                     family['elemFamilySize'],
                     family['elemFamilyType'],
                     family['elemFamilyMatId'],
                     family['ptrElemFamilyNodes'],
                     self.elemIntegrationSolid['numElemNodes'])
      
      # Get material model for family, and initialize properties using spatial database.
      #*******  Still need to figure out how to do this *********

      # Define function pointers for material-specific routines.

      # Set up state variable storage for family, based on material type and integration info.

      
    def __init__(self, name="elemfamilysolid"):
      """Constructor."""
      ElemFamily.__init__(self, name, facility="elemfamilysolid")

      self.label = ""
      self.elemType = 0
      self.quadratureOrder = ""
      self.elemInfoSolid = {'': None}
      self.stateVars = {'numStateVars': 0,
                        'numState0Vars': 0,
                        'sizeStateVars': 0,
                        'sizeState0Vars': 0,
                        'ptrStateVars': None,
                        'ptrState0Vars': None}
      # Still not quite sure how this should work.  This setup assumes
      # that I use the spatial database to put everything into a material
      # property array.
      self.matProps = {'dimProps': 0,
                       'sizeProps': 0,
                       'ptrProps': None}
      # Things that should be in ElemFamily:
      # ElemInfoSolid:  dimension and integration info
      # ElemFamilyState:stresses, strains, etc.
      # Materials:      material types (including all necessary function pointers),
      #                 database info?, property array.
      # Label:          descriptor for family
      # elemType:       string or integer description.
      # quadratureOrder:global for family

      import journal
      self._info = journal.info(name)
      return

    # PRIVATE METHODS /////////////////////////////////////////////////////////

    def _configure(self):
      """Setup members using inventory."""
      self.label = self.inventory.label
      self.quadratureOrder = self.inventory.quadratureOrder
      self.material = self.inventory.material
      self.meshimporter = self.inventory.meshimporter
      return

# version
__id__ = "$Id$"

# End of file 
