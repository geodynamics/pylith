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
 * @file unittests/libtests/feassemble/TestQuadratureEngine.hh
 *
 * @brief C++ TestQuadratureEngine object
 *
 * C++ unit testing for Quadrature.
 */

#if !defined(pylith_feassemble_testquadrature_hh)
#define pylith_feassemble_testquadrature_hh

#include <cppunit/extensions/HelperMacros.h>

#include "pylith/feassemble/feassemblefwd.hh" // forward declarations

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestQuadratureEngine;
    class QuadratureData;
  } // feassemble
} // pylith

/// C++ unit testing for Quadrature
class pylith::feassemble::TestQuadratureEngine : public CppUnit::TestFixture
{ // class TestQuadratureEngine

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestQuadratureEngine );

  CPPUNIT_TEST( testCopyConstructor );
  CPPUNIT_TEST( testCheckConditioning );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test copy constructor.
  void testCopyConstructor(void);

  /// Test checkConditioning()
  void testCheckConditioning(void);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Test computeGeometry() and retrieveGeometry().
   *
   * @param engine Quadrature engine.
   * @param data Data for testing quadrature
   */
  void _testComputeGeometry(QuadratureEngine* engine,
			    const QuadratureData& data) const;

}; // class TestQuadratureEngine

#endif // pylith_feassemble_testquadratureengine_hh

// End of file 
