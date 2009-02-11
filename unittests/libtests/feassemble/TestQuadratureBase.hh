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
 * @file unittests/libtests/feassemble/TestQuadratureBase.hh
 *
 * @brief C++ TestQuadratureBase object
 *
 * C++ unit testing for QuadratureBase.
 */

#if !defined(pylith_feassemble_testquadraturebase_hh)
#define pylith_feassemble_testquadraturebase_hh

#include <cppunit/extensions/HelperMacros.h>

/// Namespace for pylith package
namespace pylith {
  namespace feassemble {
    class TestQuadratureBase;
    class QuadratureBaseData;
  } // feassemble
} // pylith

/// C++ unit testing for QuadratureBase
class pylith::feassemble::TestQuadratureBase : public CppUnit::TestFixture
{ // class TestQuadratureBase

  // CPPUNIT TEST SUITE /////////////////////////////////////////////////
  CPPUNIT_TEST_SUITE( TestQuadratureBase );

  CPPUNIT_TEST( testConstructor );
  CPPUNIT_TEST( testMinJacobian );
  CPPUNIT_TEST( testRefGeometry );
  CPPUNIT_TEST( testInitialize );

  CPPUNIT_TEST_SUITE_END();

  // PUBLIC METHODS /////////////////////////////////////////////////////
public :

  /// Test constructor.
  void testConstructor(void);

  /// Test minJacobian()
  void testMinJacobian(void);

  /// Test refGeometry()
  void testRefGeometry(void);

  /// Test initialize()
  void testInitialize(void);

}; // class TestQuadratureBase

#endif // pylith_feassemble_testquadraturebase_hh

// End of file 
