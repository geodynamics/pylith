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

# ElemFamily class
class ElemFamily(Component):
  """Python manager for elements associated with a particular material and
  element type."""

  # INVENTORY /////////////////////////////////////////////////////////////////

  class Inventory(Component.Inventory):
    """Python object for managing ElemFamily facilities and properties."""

    import pyre.inventory

    label = pyre.inventory.str("label", default="")
    label.meta['tip'] = "Label for element family."

    # At present, it is not possible to mix different element types, but I
    # still think it makes sense to put integration info into the element
    # family.  Also, I plan to put back in the ability to mix elements.
    from Integrator import Integrator
    integrator = pyre.inventory.facility("integrator", factory=Integrator)
    integrator.meta['tip'] = "Integration information for elements."

    from Material import Material
    material = pyre.inventory.facility("material", factory=Material)
    material.meta['tip'] = "Material for a group of elements."
    
    # PUBLIC METHODS //////////////////////////////////////////////////////////

    def __init__(self, name="elemfamily"):
      """Constructor."""
      Component.__init__(self, name, facility="elemfamily")

      self.label = None
      self.elemType = None
      self.integrator = None
      self.material = None
      
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
