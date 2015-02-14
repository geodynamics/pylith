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

__all__ = ['PhysicalProperties',
           'Solution',
           'StateVariables',
           'Fault',
           ]


def has_h5py():
  if not "flag" in dir(has_h5py):
    try:
      import h5py
      has_h5py.flag = True
    except ImportError:
      print "WARNING: Cannot find h5py Python modele."
      print "         Tests limited to running PyLith without errors."
      print "         Install h5py (available via the installer utility) "
      print "         in order to enable verification of output."
      has_h5py.flag = False
  return has_h5py.flag


# End of file
