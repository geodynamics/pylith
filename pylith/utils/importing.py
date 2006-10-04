#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# {LicenseText}
#
# ----------------------------------------------------------------------
#

## @file pylith/utils/importing.py

def importModule(moduleName):
  """
  Import the named module, and return the module object
  - Works properly for fully qualified names.
  """
  module     = __import__(moduleName)
  components = moduleName.split('.')
  for comp in components[1:]:
    module = getattr(module, comp)
  return module


# End of file
