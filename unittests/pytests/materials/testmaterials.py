#!/usr/bin/env python
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

## @file unittests/materials/testmaterials.py

## @brief Python application for testing materials code.

from pyre.applications.Script import Script

import unittest

class TestApp(Script):
  """
  Test application.
  """

  # PUBLIC METHODS /////////////////////////////////////////////////////

  def __init__(self, name="testapp"):
    """
    Constructor.
    """
    Script.__init__(self, name)
    return


  def main(self):
    """
    Run the application.
    """
    from pylith.utils.PetscManager import PetscManager
    petsc = PetscManager()
    petsc.initialize()

    unittest.TextTestRunner(verbosity=2).run(self._suite())

    petsc.finalize()
    return


  # PRIVATE METHODS ////////////////////////////////////////////////////

  def _suite(self):
    """
    Setup the test suite.
    """

    suite = unittest.TestSuite()

    from TestElasticStrain1D import TestElasticStrain1D
    suite.addTest(unittest.makeSuite(TestElasticStrain1D))

    from TestElasticStress1D import TestElasticStress1D
    suite.addTest(unittest.makeSuite(TestElasticStress1D))

    from TestElasticPlaneStrain import TestElasticPlaneStrain
    suite.addTest(unittest.makeSuite(TestElasticPlaneStrain))

    from TestElasticPlaneStress import TestElasticPlaneStress
    suite.addTest(unittest.makeSuite(TestElasticPlaneStress))

    from TestElasticIsotropic3D import TestElasticIsotropic3D
    suite.addTest(unittest.makeSuite(TestElasticIsotropic3D))

    #from TestMaxwellIsotropic3D import TestMaxwellIsotropic3D
    #suite.addTest(unittest.makeSuite(TestMaxwellIsotropic3D))

    #from TestGenMaxwellIsotropic3D import TestGenMaxwellIsotropic3D
    #suite.addTest(unittest.makeSuite(TestGenMaxwellIsotropic3D))

    from TestMaterial import TestMaterial
    suite.addTest(unittest.makeSuite(TestMaterial))

    from TestElasticMaterial import TestElasticMaterial
    suite.addTest(unittest.makeSuite(TestElasticMaterial))

    from TestHomogeneous import TestHomogeneous
    suite.addTest(unittest.makeSuite(TestHomogeneous))

    return suite


# ----------------------------------------------------------------------
if __name__ == '__main__':
  app = TestApp()
  app.run()


# End of file 
