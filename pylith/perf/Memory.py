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
# Copyright (c) 2010 University of California, Davis
#
# See COPYING for license information.
#
# ----------------------------------------------------------------------
#

## @file pylith/perf/Memory.py
##
## @brief Python object for calculating memory usage.

import pylith.utils.petsc as petsc

class Memory(object):
  sizeInt    = 4
  sizeDouble = 8
  import distutils.sysconfig
  pointerSize = distutils.sysconfig.get_config_var('SIZEOF_VOID_P')
  if pointerSize is None:
    # Get pointer using sizeof(void*) in PyLith C++ library.
    pointerSize = petsc.sizeofVoidPtr()

  if pointerSize == 4:
    sizeSetEntry = 12
    sizeMapEntry = 16
    sizeArrow    = 40 # 32 bit

  elif pointerSize == 8:
    sizeSetEntry = 24
    sizeMapEntry = 32
    sizeArrow    = 56 # 64 bit

  elif pointerSize is None:
    sizeSetEntry = 0
    sizeMapEntry = 0
    sizeArrow    = 0 # Use 0 if can't get estimate of pointer size.
  else:
    raise RuntimeError("Unknown pointer size '%s'.")


# End of file
