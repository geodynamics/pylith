#!/usr/bin/env python
#
# ----------------------------------------------------------------------
#
# Brad T. Aagaard, U.S. Geological Survey
# Charles A. Williams, GNS Science
# Matthew G. Knepley, University of Chicago
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2015 University of California, Davis
#
# See COPYING for license information.
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
