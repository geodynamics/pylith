#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                             Charles A. Williams
#                       Rensselaer Polytechnic Institute
#                        (C) 2004  All Rights Reserved
#
#  Copyright 2004 Rensselaer Polytechnic Institute.
#  All worldwide rights reserved.  A license to use, copy, modify and
#  distribute this software for non-commercial research purposes only
#  is hereby granted, provided that this copyright notice and
#  accompanying disclaimer is not modified or removed from the software.
#
#  DISCLAIMER:  The software is distributed "AS IS" without any express
#  or implied warranty, including but not limited to, any implied
#  warranties of merchantability or fitness for a particular purpose
#  or any warranty of non-infringement of any current or pending patent
#  rights.  The authors of the software make no representations about
#  the suitability of this software for any particular purpose.  The
#  entire risk as to the quality and performance of the software is with
#  the user.  Should the software prove defective, the user assumes the
#  cost of all necessary servicing, repair or correction.  In
#  particular, neither Rensselaer Polytechnic Institute, nor the authors
#  of the software are liable for any indirect, special, consequential,
#  or incidental damages related to the software, to the maximum extent
#  the law permits.
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
      props.str("asciiOutputFile", default="lithomop_solver.ascii"),
      props.str("ucdOutputFile", default="lithomop_solver.inp")
      ]

# version
__id__ = "$Id: SolverOutputAsciiUcd.py,v 1.1 2004/09/23 17:37:55 willic3 Exp $"

# End of file 
