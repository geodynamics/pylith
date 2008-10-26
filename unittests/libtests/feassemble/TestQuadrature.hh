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
 * @file unittests/libtests/feassemble/TestQuadrature.hh
 *
 * @brief C++ TestQuadrature object
 *
 * C++ unit testing for Quadrature.
 */

#if !defined(pylith_feassemble_testquadrature_hh)
#define pylith_feassemble_testquadrature_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class Quadrature;
    class TestQuadrature;
    class QuadratureData;
  } // feassemble
} // pylith

/// C++ unit testing for Quadrature
class pylith::feassemble::TestQuadrature : public CppUnit::TestFixture
{ // class TestQuadrature

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestQuadrature );

  CPPUNIT_TEST( testClone );
  CPPUNIT_TEST( testMinJacobian );
  CPPUNIT_TEST( testCheckConditioning );
  CPPUNIT_TEST( testRefGeometry );
  CPPUNIT_TEST( testInitialize );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test clone()
  void testClone(void);

  /// Test minJacobian()
  void testMinJacobian(void);

  void testCheckConditioning(void);
  /// Test checkConditioning()

  /// Test refGeometry()
  void testRefGeometry(void);

  /// Test initialize()
  void testInitialize(void);

  // PROTECTED METHODS //////////////////////////////////////////////////
protected :

  /** Test initialize() & computeGeometry()
   *
   * @param pQuad Pointer to quadrature
   * @param data Data for testing quadrature
   */
  void _testComputeGeometry(Quadrature* pQuad,
			    const QuadratureData& data) const;

}; // class TestQuadrature

#endif // pylith_feassemble_testquadrature_hh

// End of file 
