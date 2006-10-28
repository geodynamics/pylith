// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
//                           Brad T. Aagaard
//                        U.S. Geological Survey
//
// {LicenseText}
//
// ----------------------------------------------------------------------
//

/**
 * @file unittests/libtests/feassemble/TestIntegrator.hh
 *
 * @brief C++ TestIntegrator object
 *
 * C++ unit testing for Integrator.
 */

#if !defined(pylith_feassemble_testintegrator_hh)
#define pylith_feassemble_testintegrator_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for spatialdata package
namespace pylith {
  namespace feassemble {
    class Integrator;
    class TestIntegrator;
    class IntegratorData;
  } // feassemble
} // pylith

/// C++ unit testing for Integrator
class pylith::feassemble::TestIntegrator : public CppUnit::TestFixture
{ // class TestIntegrator

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestIntegrator );
  CPPUNIT_TEST( testClone );
  CPPUNIT_TEST( testQuadrature );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test clone()
  void testClone(void);

  /// Test quadrature()
  void testQuadrature(void);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Test integrateAction()
   *
   * @param integrator Pointer to integrator
   * @param data Data for testing integrator
   */
  void _testIntegrateAction(Integrator* integrator,
			    const IntegratorData& data) const;

  /** Test integrate()
   *
   * @param integrator Pointer to integrator
   * @param data Data for testing integrator
   */
  void _testIntegrate(Integrator* integrator,
		      const IntegratorData& data) const;

}; // class TestIntegrator

#endif // pylith_feassemble_testintegrator_hh

// End of file 
