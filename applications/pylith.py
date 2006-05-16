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

## @file applications/pylith.py
## @brief Python PyLith application driver

from pylith.problems.PyLithApp import PyLithApp

if __name__ == "__main__":
  app = PyLithApp()
  app.run()


# version
__id__ = "$Id$"

# End of file 
