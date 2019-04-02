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

from pylith.utils.utils import PetscVersion

class TestPetscVersion(unittest.TestCase):

  def test_isRelease(self):
    isRelease = PetscVersion.isRelease()
    return


  def test_version(self):
    version = PetscVersion.version()
    # Check that version is of the form X.X.X
    import re
    match = re.search("[0-9]+\.[0-9]+\.[0-9]+", version)
    self.failIf(match is None)
    return


  def test_gitVersion(self):
    revision = PetscVersion.gitRevision()
    if PetscVersion.isRelease():
      pass
    else:
      # Check that revision is of the form v2.1.3-16-g9323114 or v3.10-88-g06a760874e
      import re
      match = re.search("v[0-9]+\.[0-9]+(\.[0-9]+)*-[0-9]+-g[0-9,a-z]+", revision)
      self.failIf(match is None)
    return


  def test_gitDate(self):
    value = PetscVersion.gitDate()
    if PetscVersion.isRelease():
      pass
    else:
      # Check form of datetime
      import datetime
      fields = value.split()
      d = datetime.datetime.strptime(fields[0], "%Y-%m-%d")
      t = datetime.datetime.strptime(fields[1], "%H:%M:%S")
    return


  def test_gitBranch(self):
    branch = PetscVersion.gitBranch()
    if PetscVersion.isRelease():
      pass
    else:
      self.failIf(len(branch) == 0)
    return


  def test_petscDir(self):
    dir = PetscVersion.petscDir()
    self.failIf(len(dir) == 0)
    return


  def test_petscArch(self):
    arch = PetscVersion.petscArch()
    # Prefix builds set PETSC_ARCH to empty string, so no verification of string length.
    return


# End of file 
