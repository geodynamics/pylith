#!/usr/bin/env python
#
# ======================================================================
#
# Brad T. Aagaard, U.S. Geological Survey
#
# This code was developed as part of the Computational Infrastructure
# for Geodynamics (http://geodynamics.org).
#
# Copyright (c) 2010-2017 University of California, Davis
#
# See COPYING for license information.
#
# ======================================================================
#

import unittest

from pylith.utils.utils import DependenciesVersion

class TestDependenciesVersion(unittest.TestCase):

  def test_mpiVersion(self):
    version = DependenciesVersion.mpiVersion()
    # Check that version is of the form X.X.X or X.X
    import re
    match = re.search("[0-9]+\.[0-9]+\.[0-9]+", version)
    if match is None:
      match = re.search("[0-9]+\.[0-9]+", version)
    self.failIf(match is None)
    return


  def test_mpiImplementation(self):
    imp = DependenciesVersion.mpiImplementation()
    self.failIf(len(imp) == 0)
    return


  def test_mpiStandard(self):
    version = DependenciesVersion.mpiStandard()
    # Check that version is of the form X.X
    import re
    match = re.search("[0-9]+\.[0-9]+", version)
    self.failIf(match is None)
    return


  def test_netcdfVersion(self):
    version = DependenciesVersion.netcdfVersion()
    # Check that version is of the form X.X.X
    import re
    match = re.search("[0-9]+\.[0-9]+\.[0-9]+", version)
    self.failIf(match is None)
    return


  def test_hdf5Version(self):
    version = DependenciesVersion.hdf5Version()
    # Check that version is of the form X.X.X
    import re
    match = re.search("[0-9]+\.[0-9]+\.[0-9]+", version)
    self.failIf(match is None)
    return


# End of file 
