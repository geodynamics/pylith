# -*- Pyrex -*-
#
# ======================================================================
#
#                           Brad T. Aagaard
#                        U.S. Geological Survey
#
# {LicenseText}
#
# ======================================================================
#

# ----------------------------------------------------------------------
cdef class Mesh:

  cdef void* thisptr # Pointer to C++ object
  cdef readonly object handle # PyCObject holding pointer to C++ object
  cdef readonly object name # Identifier for object base type


# End of file 
