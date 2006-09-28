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

from SolverOutput import SolverOutput

# ----------------------------------------------------------------------
class SolverOutputAsciiUcd(SolverOutput):

  def print_status(self):

    status = journal.info("Output.info")
    status.line("  filename: %s" % self.inventory.filename)
    status.log()
    return

  def _writeTimeStep(self, t, XXXXXXXX):

    if self.inventory.asciiOutput:
      bindings.solver_output_ascii(XXXXXX)
    if self.inventory.ucdOutput:
      bindings.solver_output_ucd(XXXXXX)
    return

  def __init__(self):

    from pyre.units.time import second

    Component.__init__(self, "OutputAsciiUcd", "OutputAsciiUcd")

    self._progress = journal.debug("Output.progress")
    return

  class Inventory(Component.Inventory):
      
    import pyre.properties as props
    from pyre.units.time import year
    
    inventory = [
      props.bool("asciiOutput", default = False),
      props.bool("ucdOutput", default=True),
      props.str("asciiOutputFile", default="pylith_solver.ascii"),
      props.str("ucdOutputFile", default="pylith_solver.inp")
      ]

# version
__id__ = "$Id: SolverOutputAsciiUcd.py,v 1.1 2004/09/23 17:37:55 willic3 Exp $"

# End of file 
