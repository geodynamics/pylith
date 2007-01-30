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

from pylith.PyLithApp import PyLithApp

if __name__ == "__main__":

  import journal
  journal.info('coordinates').activate()

  app = PyLithApp()
  app.run()


# End of file 
